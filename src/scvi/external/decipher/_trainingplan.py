import pyro
import torch

from scvi.module.base import (
    PyroBaseModuleClass,
)
from scvi.train import LowLevelPyroTrainingPlan


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
            optim_kwargs.update({"lr": 5e-3})
        if "weight_decay" not in optim_kwargs.keys():
            optim_kwargs.update({"weight_decay": 1e-4})
        self.optim = pyro.optim.ClippedAdam(optim_args=optim_kwargs) if optim is None else optim
        # We let SVI take care of all optimization
        self.automatic_optimization = False

        self.svi = pyro.infer.SVI(
            model=self.module.model,
            guide=self.module.guide,
            optim=self.optim,
            loss=self.loss_fn,
        )
        self.validation_step_outputs = []

        # See configure_optimizers for what this does
        self._dummy_param = torch.nn.Parameter(torch.Tensor([0.0]))

    def on_validation_model_train(self):
        """Prepare the model for validation by switching to train mode."""
        super().on_validation_model_train()
        if self.current_epoch > 0:
            # freeze the batch norm layers after the first epoch
            # 1) the batch norm layers helps with the initialization
            # 2) but then, they seem to imply a strong normal prior on the latent space
            for module in self.module.modules():
                if isinstance(module, torch.nn.BatchNorm1d):
                    module.eval()

    def training_step(self, batch, batch_idx):
        """Training step for Pyro training."""
        args, kwargs = self.module._get_fn_args_from_batch(batch)

        # pytorch lightning requires a Tensor object for loss
        loss = torch.Tensor([self.svi.step(*args, **kwargs)])
        n_obs = args[0].shape[0]

        _opt = self.optimizers()
        _opt.step()

        out_dict = {"loss": loss, "n_obs": n_obs}
        self.training_step_outputs.append(out_dict)
        return out_dict

    def on_train_epoch_end(self):
        """Training epoch end for Pyro training."""
        outputs = self.training_step_outputs
        elbo = 0
        n_obs = 0
        for out in outputs:
            elbo += out["loss"]
            n_obs += out["n_obs"]
        elbo /= n_obs
        self.log("elbo_train", elbo, prog_bar=True)
        self.training_step_outputs.clear()

    def validation_step(self, batch, batch_idx):
        """Validation step for Pyro training."""
        args, kwargs = self.module._get_fn_args_from_batch(batch)
        loss = self.differentiable_loss_fn(
            self.scale_fn(self.module.model),
            self.scale_fn(self.module.guide),
            *args,
            **kwargs,
        )
        nll = -self.module.predictive_log_likelihood(*args, **kwargs, n_samples=5)
        out_dict = {"loss": loss, "nll": nll, "n_obs": args[0].shape[0]}
        self.validation_step_outputs.append(out_dict)
        return out_dict

    def on_validation_epoch_end(self):
        """Validation epoch end for Pyro training."""
        outputs = self.validation_step_outputs
        elbo = 0
        nll = 0
        n_obs = 0
        for out in outputs:
            elbo += out["loss"]
            nll += out["nll"]
            n_obs += out["n_obs"]
        elbo /= n_obs
        nll /= n_obs
        self.log("elbo_validation", elbo, prog_bar=True)
        self.log("nll_validation", nll, prog_bar=False)
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
