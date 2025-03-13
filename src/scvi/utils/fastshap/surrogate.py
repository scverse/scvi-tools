from copy import deepcopy

import numpy as np
import torch
import torch.optim as optim
from torch.utils.data import BatchSampler, DataLoader, Dataset, RandomSampler, TensorDataset
from tqdm.auto import tqdm

from scvi.utils.fastshap.utils import DatasetRepeat, UniformSampler


def validate(surrogate, loss_fn, data_loader):
    """
    Calculate mean validation loss.

    Args:
      loss_fn: loss function.
      data_loader: data loader.
    """
    with torch.no_grad():
        # Setup.
        device = next(surrogate.surrogate.parameters()).device
        mean_loss = 0
        N = 0

        for x, y, S in data_loader:
            x = x.to(device)
            y = y.to(device)
            S = S.to(device)
            pred = surrogate(x, S)
            loss = loss_fn(pred, y)
            N += len(x)
            mean_loss += len(x) * (loss - mean_loss) / N

    return mean_loss


def generate_labels(dataset, model, batch_size):
    """
    Generate prediction labels for a set of inputs.

    Args:
      dataset: dataset object.
      model: predictive model.
      batch_size: minibatch size.
    """
    with torch.no_grad():
        # Setup.
        preds = []
        if isinstance(model, torch.nn.Module):
            device = next(model.parameters()).device
        else:
            device = torch.device("cpu")
        loader = DataLoader(dataset, batch_size=batch_size)

        for (x,) in loader:
            pred = model(x.to(device))
            if device.type != "cpu":
                pred = pred.cpu()
            preds.append(pred)

    if type(preds[0]).__name__ == "DataFrame":
        preds = [torch.tensor(preds[0].values, dtype=torch.float32, device=device)]

    return torch.cat(preds)


class Surrogate:
    """
    Wrapper around surrogate model.

    Args:
      surrogate: surrogate model.
      num_features: number of features.
      groups: (optional) feature groups, represented by a list of lists.
    """

    def __init__(self, surrogate, num_features, groups=None):
        # Store surrogate model.
        self.surrogate = surrogate

        # Store feature groups.
        if groups is None:
            self.num_players = num_features
            self.groups_matrix = None
        else:
            # Verify groups.
            inds_list = []
            for group in groups:
                inds_list += list(group)
            assert np.all(np.sort(inds_list) == np.arange(num_features))

            # Map groups to features.
            self.num_players = len(groups)
            device = next(surrogate.parameters()).device
            self.groups_matrix = torch.zeros(
                len(groups), num_features, dtype=torch.float32, device=device
            )
            for i, group in enumerate(groups):
                self.groups_matrix[i, group] = 1

    def train(
        self,
        train_data,
        val_data,
        batch_size,
        max_epochs,
        loss_fn,
        validation_samples=1,
        validation_batch_size=None,
        lr=1e-3,
        min_lr=1e-5,
        lr_factor=0.5,
        lookback=5,
        training_seed=None,
        validation_seed=None,
        bar=False,
        verbose=False,
    ):
        """
        Train surrogate model.

        Args:
          train_data: training data with inputs and the original model's
            predictions (np.ndarray tuple, torch.Tensor tuple,
            torch.utils.data.Dataset).
          val_data: validation data with inputs and the original model's
            predictions (np.ndarray tuple, torch.Tensor tuple,
            torch.utils.data.Dataset).
          batch_size: minibatch size.
          max_epochs: maximum training epochs.
          loss_fn: loss function (e.g., fastshap.KLDivLoss).
          validation_samples: number of samples per validation example.
          validation_batch_size: validation minibatch size.
          lr: initial learning rate.
          min_lr: minimum learning rate.
          lr_factor: learning rate decrease factor.
          lookback: lookback window for early stopping.
          training_seed: random seed for training.
          validation_seed: random seed for generating validation data.
          verbose: verbosity.
        """
        # Set up train dataset.
        if isinstance(train_data, tuple):
            x_train, y_train = train_data
            if isinstance(x_train, np.ndarray):
                x_train = torch.tensor(x_train, dtype=torch.float32)
                y_train = torch.tensor(y_train, dtype=torch.float32)
            train_set = TensorDataset(x_train, y_train)
        elif isinstance(train_data, Dataset):
            train_set = train_data
        else:
            raise ValueError("train_data must be either tuple of tensors or a PyTorch Dataset")

        # Set up train data loader.
        random_sampler = RandomSampler(
            train_set,
            replacement=True,
            num_samples=int(np.ceil(len(train_set) / batch_size)) * batch_size,
        )
        batch_sampler = BatchSampler(random_sampler, batch_size=batch_size, drop_last=True)
        train_loader = DataLoader(train_set, batch_sampler=batch_sampler)

        # Set up validation dataset.
        sampler = UniformSampler(self.num_players)
        if validation_seed is not None:
            torch.manual_seed(validation_seed)
        S_val = sampler.sample(len(val_data) * validation_samples)

        if isinstance(val_data, tuple):
            x_val, y_val = val_data
            if isinstance(x_val, np.ndarray):
                x_val = torch.tensor(x_val, dtype=torch.float32)
                y_val = torch.tensor(y_val, dtype=torch.float32)
            x_val_repeat = x_val.repeat(validation_samples, 1)
            y_val_repeat = y_val.repeat(validation_samples, 1)
            val_set = TensorDataset(x_val_repeat, y_val_repeat, S_val)
        elif isinstance(val_data, Dataset):
            val_set = DatasetRepeat([val_data, TensorDataset(S_val)])
        else:
            raise ValueError("val_data must be either tuple of tensors or a PyTorch Dataset")

        if validation_batch_size is None:
            validation_batch_size = batch_size
        val_loader = DataLoader(val_set, batch_size=validation_batch_size)

        # Setup for training.
        surrogate = self.surrogate
        device = next(surrogate.parameters()).device
        optimizer = optim.Adam(surrogate.parameters(), lr=lr)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, factor=lr_factor, patience=lookback // 2, min_lr=min_lr, verbose=verbose
        )
        best_loss = validate(self, loss_fn, val_loader).item()
        best_epoch = 0
        best_model = deepcopy(surrogate)
        loss_list = [best_loss]
        if training_seed is not None:
            torch.manual_seed(training_seed)

        for epoch in range(max_epochs):
            # Batch iterable.
            if bar:
                batch_iter = tqdm(train_loader, desc="Training epoch")
            else:
                batch_iter = train_loader

            for x, y in batch_iter:
                # Prepare data.
                x = x.to(device)
                y = y.to(device)

                # Generate subsets.
                S = sampler.sample(batch_size).to(device=device)

                # Make predictions.
                pred = self.__call__(x, S)
                loss = loss_fn(pred, y)

                # Optimizer step.
                loss.backward()
                optimizer.step()
                surrogate.zero_grad()

            # Evaluate validation loss.
            self.surrogate.eval()
            val_loss = validate(self, loss_fn, val_loader).item()
            self.surrogate.train()

            # Print progress.
            if verbose:
                print(f"----- Epoch = {epoch + 1} -----")
                print(f"Val loss = {val_loss:.4f}")
                print("")
            scheduler.step(val_loss)
            loss_list.append(val_loss)

            # Check if best model.
            if val_loss < best_loss:
                best_loss = val_loss
                best_model = deepcopy(surrogate)
                best_epoch = epoch
                if verbose:
                    print(f"New best epoch, loss = {val_loss:.4f}")
                    print("")
            elif epoch - best_epoch == lookback:
                if verbose:
                    print("Stopping early")
                break

        # Clean up.
        for param, best_param in zip(
            surrogate.parameters(), best_model.parameters(), strict=False
        ):
            param.data = best_param.data
        self.loss_list = loss_list
        self.surrogate.eval()

    def train_original_model(
        self,
        train_data,
        val_data,
        original_model,
        batch_size,
        max_epochs,
        loss_fn,
        validation_samples=1,
        validation_batch_size=None,
        lr=1e-3,
        min_lr=1e-5,
        lr_factor=0.5,
        lookback=5,
        training_seed=None,
        validation_seed=None,
        bar=False,
        verbose=False,
    ):
        """
        Train surrogate model with labels provided by the original model.

        Args:
          train_data: training data with inputs only (np.ndarray, torch.Tensor,
            torch.utils.data.Dataset).
          val_data: validation data with inputs only (np.ndarray, torch.Tensor,
            torch.utils.data.Dataset).
          original_model: original predictive model (e.g., torch.nn.Module).
          batch_size: minibatch size.
          max_epochs: maximum training epochs.
          loss_fn: loss function (e.g., fastshap.KLDivLoss).
          validation_samples: number of samples per validation example.
          validation_batch_size: validation minibatch size.
          lr: initial learning rate.
          min_lr: minimum learning rate.
          lr_factor: learning rate decrease factor.
          lookback: lookback window for early stopping.
          training_seed: random seed for training.
          validation_seed: random seed for generating validation data.
          verbose: verbosity.
        """
        # Set up train dataset.
        if isinstance(train_data, np.ndarray):
            train_data = torch.tensor(train_data, dtype=torch.float32)

        if isinstance(train_data, torch.Tensor):
            train_set = TensorDataset(train_data)
        elif isinstance(train_data, Dataset):
            train_set = train_data
        else:
            raise ValueError("train_data must be either tensor or a PyTorch Dataset")

        # Set up train data loader.
        random_sampler = RandomSampler(
            train_set,
            replacement=True,
            num_samples=int(np.ceil(len(train_set) / batch_size)) * batch_size,
        )
        batch_sampler = BatchSampler(random_sampler, batch_size=batch_size, drop_last=True)
        train_loader = DataLoader(train_set, batch_sampler=batch_sampler)

        # Set up validation dataset.
        sampler = UniformSampler(self.num_players)
        if validation_seed is not None:
            torch.manual_seed(validation_seed)
        S_val = sampler.sample(len(val_data) * validation_samples)
        if validation_batch_size is None:
            validation_batch_size = batch_size

        if isinstance(val_data, np.ndarray):
            val_data = torch.tensor(val_data, dtype=torch.float32)

        if isinstance(val_data, torch.Tensor):
            # Generate validation labels.
            y_val = generate_labels(TensorDataset(val_data), original_model, validation_batch_size)
            y_val_repeat = y_val.repeat(validation_samples, *[1 for _ in y_val.shape[1:]])

            # Create dataset.
            val_data_repeat = val_data.repeat(validation_samples, 1)
            val_set = TensorDataset(val_data_repeat, y_val_repeat, S_val)
        elif isinstance(val_data, Dataset):
            # Generate validation labels.
            y_val = generate_labels(val_data, original_model, validation_batch_size)
            y_val_repeat = y_val.repeat(validation_samples, *[1 for _ in y_val.shape[1:]])

            # Create dataset.
            val_set = DatasetRepeat([val_data, TensorDataset(y_val_repeat, S_val)])
        else:
            raise ValueError("val_data must be either tuple of tensors or a PyTorch Dataset")

        val_loader = DataLoader(val_set, batch_size=validation_batch_size)

        # Setup for training.
        surrogate = self.surrogate
        device = next(surrogate.parameters()).device
        optimizer = optim.Adam(surrogate.parameters(), lr=lr)
        scheduler = optim.lr_scheduler.ReduceLROnPlateau(
            optimizer, factor=lr_factor, patience=lookback // 2, min_lr=min_lr, verbose=verbose
        )
        best_loss = validate(self, loss_fn, val_loader).item()
        best_epoch = 0
        best_model = deepcopy(surrogate)
        loss_list = [best_loss]
        if training_seed is not None:
            torch.manual_seed(training_seed)

        for epoch in range(max_epochs):
            # Batch iterable.
            if bar:
                batch_iter = tqdm(train_loader, desc="Training epoch")
            else:
                batch_iter = train_loader

            for (x,) in batch_iter:
                # Prepare data.
                x = x.to(device)

                # Get original model prediction.
                with torch.no_grad():
                    y = original_model(x.cpu())

                # Generate subsets.
                S = sampler.sample(batch_size).to(device=device)

                # Make predictions.
                pred = self.__call__(x, S)
                if type(y).__name__ == "DataFrame":
                    y = torch.tensor(y.values, dtype=torch.float32, device=device)
                loss = loss_fn(pred, y)

                # Optimizer step.
                loss.backward()
                optimizer.step()
                surrogate.zero_grad()

            # Evaluate validation loss.
            self.surrogate.eval()
            val_loss = validate(self, loss_fn, val_loader).item()
            self.surrogate.train()

            # Print progress.
            if verbose:
                print(f"----- Epoch = {epoch + 1} -----")
                print(f"Val loss = {val_loss:.4f}")
                print("")
            scheduler.step(val_loss)
            loss_list.append(val_loss)

            # Check if best model.
            if val_loss < best_loss:
                best_loss = val_loss
                best_model = deepcopy(surrogate)
                best_epoch = epoch
                if verbose:
                    print(f"New best epoch, loss = {val_loss:.4f}")
                    print("")
            elif epoch - best_epoch == lookback:
                if verbose:
                    print("Stopping early")
                break

        # Clean up.
        for param, best_param in zip(
            surrogate.parameters(), best_model.parameters(), strict=False
        ):
            param.data = best_param.data
        self.loss_list = loss_list
        self.surrogate.eval()

    def __call__(self, x, S):
        """
        Evaluate surrogate model.

        Args:
          x: input examples.
          S: coalitions.
        """
        if self.groups_matrix is not None:
            S = torch.mm(S, self.groups_matrix)

        return self.surrogate((x, S))
