import logging

import numpy as np
import torch
from sklearn.neighbors import KNeighborsClassifier

from scvi import _CONSTANTS
from scvi.core import unsupervised_clustering_accuracy

from .scvi_data_loader import ScviDataLoader
from scvi.core._log_likelihood import compute_elbo

logger = logging.getLogger(__name__)


class AnnotationDataLoader(ScviDataLoader):
    def __init__(self, *args, unlabeled=False, model_zl=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.model_zl = model_zl
        self.unlabeled = unlabeled

    def accuracy(self):
        model, cls = (
            (self.sampling_model, self.model)
            if hasattr(self, "sampling_model")
            else (self.model, None)
        )
        acc = compute_accuracy(model, self, classifier=cls, model_zl=self.model_zl)
        logger.debug("Acc: %.4f" % (acc))
        return acc

    accuracy.mode = "max"

    @torch.no_grad()
    def hierarchical_accuracy(self):
        all_y, all_y_pred = self.compute_predictions()
        acc = np.mean(all_y == all_y_pred)

        all_y_groups = np.array([self.model.labels_groups[y] for y in all_y])
        all_y_pred_groups = np.array([self.model.labels_groups[y] for y in all_y_pred])
        h_acc = np.mean(all_y_groups == all_y_pred_groups)

        logger.debug("Hierarchical Acc : %.4f\n" % h_acc)
        return acc

    accuracy.mode = "max"

    @torch.no_grad()
    def compute_predictions(self, soft=False):
        """

        Parameters
        ----------
        soft
             (Default value = False)

        Returns
        -------
        the true labels and the predicted labels

        """
        model, cls = (
            (self.sampling_model, self.model)
            if hasattr(self, "sampling_model")
            else (self.model, None)
        )
        return compute_predictions(
            model, self, classifier=cls, soft=soft, model_zl=self.model_zl
        )

    @torch.no_grad()
    def unsupervised_classification_accuracy(self):
        all_y, all_y_pred = self.compute_predictions()
        uca = unsupervised_clustering_accuracy(all_y, all_y_pred)[0]
        logger.debug("UCA : %.4f" % (uca))
        return uca

    unsupervised_classification_accuracy.mode = "max"

    @torch.no_grad()
    def elbo(self) -> torch.Tensor:
        """Returns the Evidence Lower Bound associated to the object."""
        elbo = compute_elbo(self.model, self, feed_labels=not self.unlabeled)
        logger.debug("ELBO : %.4f" % elbo)
        return elbo

    elbo.mode = "min"

    @torch.no_grad()
    def nn_latentspace(self, data_loader):
        data_train, _, labels_train = self.get_latent()
        data_test, _, labels_test = data_loader.get_latent()
        nn = KNeighborsClassifier()
        nn.fit(data_train, labels_train)
        score = nn.score(data_test, labels_test)
        return score


@torch.no_grad()
def compute_accuracy(vae, data_loader, classifier=None, model_zl=False):
    all_y, all_y_pred = compute_predictions(
        vae, data_loader, classifier=classifier, model_zl=model_zl
    )
    return np.mean(all_y == all_y_pred)


@torch.no_grad()
def compute_predictions(
    model, data_loader, classifier=None, soft=False, model_zl=False
):
    all_y_pred = []
    all_y = []

    for _, tensors in enumerate(data_loader):
        sample_batch = tensors[_CONSTANTS.X_KEY]
        batch_index = tensors[_CONSTANTS.BATCH_KEY]
        labels = tensors[_CONSTANTS.LABELS_KEY]

        all_y += [labels.view(-1).cpu()]

        if hasattr(model, "classify"):
            y_pred = model.classify(sample_batch, batch_index)
        elif classifier is not None:
            # Then we use the specified classifier
            if model is not None:
                if model.log_variational:
                    sample_batch = torch.log(1 + sample_batch)
                if model_zl:
                    sample_z = model.z_encoder(sample_batch, batch_index)[0]
                    sample_l = model.l_encoder(sample_batch, batch_index)[0]
                    sample_batch = torch.cat((sample_z, sample_l), dim=-1)
                else:
                    sample_batch, _, _ = model.z_encoder(sample_batch, batch_index)
            y_pred = classifier(sample_batch)
        else:  # The model is the raw classifier
            y_pred = model(sample_batch)

        if not soft:
            y_pred = y_pred.argmax(dim=-1)

        all_y_pred += [y_pred.cpu()]

    all_y_pred = np.array(torch.cat(all_y_pred))
    all_y = np.array(torch.cat(all_y))

    return all_y, all_y_pred
