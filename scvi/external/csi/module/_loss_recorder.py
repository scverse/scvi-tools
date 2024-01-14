class LossRecorder:
    """
    Loss signature for models.

    This class provides an organized way to record the model loss, as well as
    the components of the ELBO. This may also be used in MLE, MAP, EM methods.
    The loss is used for backpropagation during inference. The other parameters
    are used for logging/early stopping during inference.

    Parameters
    ----------
    loss
        Tensor with loss for minibatch. Should be one dimensional with one value.
        Note that loss should be a :class:`~torch.Tensor` and not the result of ``.item()``.
    reconstruction_loss
        Reconstruction loss for each observation in the minibatch.
    kl_local
        KL divergence associated with each observation in the minibatch.
    kl_global
        Global kl divergence term. Should be one dimensional with one value.
    **kwargs
        Additional metrics can be passed as keyword arguments and will
        be available as attributes of the object.
    """

    def __init__(
        self,
        n_obs: int,
        loss: float,
        loss_sum: float,
        **kwargs,
    ):
        self.n_obs = n_obs
        self.loss = loss
        self.loss_sum = loss_sum
        self.extra_metric_attrs = []
        for key, value in kwargs.items():
            setattr(self, key, value)
            self.extra_metric_attrs.append(key)
