import pyro
import torch

from scvi.module.base import (
    PyroBaseModuleClass,
)
from scvi.train import LowLevelPyroTrainingPlan

from ._utils import predictive_log_likelihood


class DecipherTrainingPlan(LowLevelPyroTrainingPlan):
    """Lightning module task to train the Decipher Pyro module.

    Parameters
    ----------
    pyro_module
        An instance of :class:`~scvi.module.base.PyroBaseModuleClass`. This object
        should have callable `model` and `guide` attributes or methods.
    loss_fn
        A Pyro loss. Should be a subclass of :class:`~pyro.infer.ELBO`.
        If `None`, defaults to :class:`~pyro.infer.Trace_ELBO`.
    optim
        A Pyro optimizer instance, e.g., :class:`~pyro.optim.Adam`. If `None`,
        defaults to :class:`pyro.optim.Adam` optimizer with a learning rate of `1e-3`.
    optim_kwargs
        Keyword arguments for **default** optimiser :class:`pyro.optim.Adam`.
    """

    def __init__(
        self,
        pyro_module: PyroBaseModuleClass,
        loss_fn: pyro.infer.ELBO | None = None,
        optim: pyro.optim.PyroOptim | None = None,
        optim_kwargs: dict | None = None,
    ):
        super().__init__(
            pyro_module=pyro_module,
            loss_fn=loss_fn,
        )
        optim_kwargs = optim_kwargs if isinstance(optim_kwargs, dict) else {}
        if "lr" not in optim_kwargs.keys():
            optim_kwargs.update({"lr": 5e-3, "weight_decay": 1e-4})
        self.optim = pyro.optim.ClippedAdam(optim_args=optim_kwargs) if optim is None else optim
        # We let SVI take care of all optimization
        self.automatic_optimization = False

        self.svi = pyro.infer.SVI(
            model=self.module.model,
            guide=self.module.guide,
            optim=self.optim,
            loss=self.loss_fn,
        )
        # See configure_optimizers for what this does
        self._dummy_param = torch.nn.Parameter(torch.Tensor([0.0]))

    def training_step(self, batch, batch_idx):
        """Training step for Pyro training."""
        args, kwargs = self.module._get_fn_args_from_batch(batch)

        # pytorch lightning requires a Tensor object for loss
        loss = torch.Tensor([self.svi.step(*args, **kwargs)])

        _opt = self.optimizers()
        _opt.step()

        out_dict = {"loss": loss}
        self.training_step_outputs.append(out_dict)
        return out_dict

    def on_train_epoch_start(self):
        """Training epoch start for Pyro training."""
        super().on_train_epoch_start()
        if self.current_epoch > 0:
            # freeze the batch norm layers after the first epoch
            # 1) the batch norm layers helps with the initialization
            # 2) but then, they seem to imply a strong normal prior on the latent space
            for module in self.module.modules():
                if isinstance(module, torch.nn.BatchNorm1d):
                    module.eval()

    def on_train_epoch_end(self):
        """Training epoch end for Pyro training."""
        outputs = self.training_step_outputs
        elbo = 0
        for out in outputs:
            elbo += out["loss"]
        elbo /= self.n_obs_training
        self.log("elbo_train", elbo, prog_bar=True)
        self.training_step_outputs.clear()

    def validation_step(self, batch, batch_idx):
        """Validation step for Pyro training."""
        out_dict = super().validation_step(batch, batch_idx)
        args, kwargs = self.module._get_fn_args_from_batch(batch)
        nll = -predictive_log_likelihood(self.module, *args, **kwargs, n_samples=5)
        out_dict["nll"] = nll
        self.validation_step_outputs[-1].update(out_dict)
        return out_dict

    def on_validation_epoch_end(self):
        """Validation epoch end for Pyro training."""
        outputs = self.validation_step_outputs
        elbo = 0
        nll = 0
        for out in outputs:
            elbo += out["loss"]
            nll += out["nll"]
        elbo /= self.n_obs_validation
        nll /= self.n_obs_validation
        self.log("elbo_validation", elbo, prog_bar=True)
        self.log("nll_validation", nll, prog_bar=True)
        self.validation_step_outputs.clear()

    def configure_optimizers(self):
        """Shim optimizer for PyTorch Lightning.

        PyTorch Lightning wants to take steps on an optimizer
        returned by this function in order to increment the global
        step count. See PyTorch Lighinting optimizer manual loop.

        Here we provide a shim optimizer that we can take steps on
        at minimal computational cost in order to keep Lightning happy :).
        """
        return torch.optim.Adam([self._dummy_param])

    def optimizer_step(self, *args, **kwargs):
        pass

    def backward(self, *args, **kwargs):
        pass
