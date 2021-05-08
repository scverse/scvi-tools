from contextlib import ExitStack  # python 3

import numpy as np
import pyro
import pyro.distributions as dist
import torch
from pyro.distributions import constraints
from pyro.distributions.transforms import SoftplusTransform
from pyro.distributions.util import sum_rightmost
from pyro.infer.autoguide import AutoGuide
from pyro.infer.autoguide import AutoGuideList as PyroAutoGuideList
from pyro.infer.autoguide.guides import _deep_getattr, _deep_setattr
from pyro.infer.autoguide.utils import helpful_support_errors
from pyro.nn import PyroModule, PyroParam
from torch.distributions import biject_to, transform_to

from scvi.nn import FCLayers


class FCLayersPyro(FCLayers, PyroModule):
    pass


class AutoGuideList(PyroAutoGuideList):
    def quantiles(self, quantiles, *args, **kwargs):
        """
        Returns the posterior quantile values of each latent variable.

        Parameters
        ----------
        quantiles
            A list of requested quantiles between 0 and 1.

        Returns
        -------
        A dict mapping sample site name to quantiles tensor.
        """

        result = {}
        for part in self:
            result.update(part.quantiles(quantiles, *args, **kwargs))
        return result


@biject_to.register(constraints.positive)
@transform_to.register(constraints.positive)
def _transform_to_positive(constraint):
    return SoftplusTransform()


class AutoNormalEncoder(AutoGuide):
    def __init__(
        self,
        model,
        amortised_plate_sites: dict,
        n_in: int,
        n_hidden: int = 200,
        init_param=0,
        init_param_scale: float = 1 / 50,
        data_transform=torch.log1p,
        encoder_class=FCLayersPyro,
        encoder_kwargs=None,
        create_plates=None,
        single_encoder: bool = True,
    ):

        super().__init__(model, create_plates=create_plates)
        self.amortised_plate_sites = amortised_plate_sites
        self.single_encoder = single_encoder

        self.softplus = SoftplusTransform()

        encoder_kwargs = encoder_kwargs if isinstance(encoder_kwargs, dict) else dict()
        encoder_kwargs["n_hidden"] = n_hidden
        self.encoder_kwargs = encoder_kwargs

        self.n_in = n_in
        self.n_out = (
            np.sum(
                [
                    np.sum(amortised_plate_sites["sites"][k])
                    for k in amortised_plate_sites["sites"].keys()
                ]
            )
            * 2
        )
        self.n_hidden = n_hidden
        self.encoder_class = encoder_class
        if self.single_encoder:
            # create a single encoder NN
            self.encoder = encoder_class(
                n_in=self.n_in, n_out=self.n_hidden, **self.encoder_kwargs
            )

        self.init_param_scale = init_param_scale
        self.data_transform = data_transform

    def _setup_prototype(self, *args, **kwargs):

        super()._setup_prototype(*args, **kwargs)

        self._event_dims = {}
        self._cond_indep_stacks = {}
        self.hidden2locs = PyroModule()
        self.hidden2scales = PyroModule()

        if not self.single_encoder:
            # create module class for collecting multiple encoder NN
            self.encoder = PyroModule()

        # Initialize guide params
        for name, site in self.prototype_trace.iter_stochastic_nodes():
            # Collect unconstrained event_dims, which may differ from constrained event_dims.
            with helpful_support_errors(site):
                init_loc = (
                    biject_to(site["fn"].support).inv(site["value"].detach()).detach()
                )
            event_dim = site["fn"].event_dim + init_loc.dim() - site["value"].dim()
            self._event_dims[name] = event_dim

            # Collect independence contexts.
            self._cond_indep_stacks[name] = site["cond_indep_stack"]

            # add linear layer for locs and scales
            param_dim = (self.n_hidden, self.amortised_plate_sites["sites"][name])
            # init_param = torch.normal(torch.full(size=param_dim, fill_value=0.),
            #                          torch.full(param_dim, 1. / np.sqrt(self.n_hidden)),
            #                          device=site["value"].device)
            init_param = np.random.normal(
                np.zeros(param_dim),
                (np.ones(param_dim) * self.init_param_scale) / np.sqrt(self.n_hidden),
            ).astype("float32")
            _deep_setattr(
                self.hidden2locs,
                name,
                PyroParam(
                    torch.tensor(
                        init_param, device=site["value"].device, requires_grad=True
                    )
                ),
            )

            init_param = np.random.normal(
                np.zeros(param_dim),
                (np.ones(param_dim) * self.init_param_scale) / np.sqrt(self.n_hidden),
            ).astype("float32")
            _deep_setattr(
                self.hidden2scales,
                name,
                PyroParam(
                    torch.tensor(
                        init_param, device=site["value"].device, requires_grad=True
                    )
                ),
            )

            if not self.single_encoder:
                _deep_setattr(
                    self.encoder,
                    name,
                    self.encoder_class(
                        n_in=self.n_in, n_out=self.n_hidden, **self.encoder_kwargs
                    ).to(site["value"].device),
                )

    def _get_loc_and_scale(self, name, encoded_hidden):

        linear_locs = _deep_getattr(self.hidden2locs, name)
        linear_scales = _deep_getattr(self.hidden2scales, name)

        if not self.single_encoder:
            # when using multiple encoders extract hidden layer for this parameter
            encoded_hidden = encoded_hidden[name]

        locs = encoded_hidden @ linear_locs
        scales = self.softplus((encoded_hidden @ linear_scales) - 2)

        return locs, scales

    def encode(self, *args, **kwargs):
        in_names = self.amortised_plate_sites["in"]
        x_in = [kwargs[i] if i in kwargs.keys() else args[i] for i in in_names]
        # apply data_transform
        x_in = [self.data_transform(x) for x in x_in]
        # when there are multiple encoders get then and apply encode data
        if not self.single_encoder:
            res = {
                name: _deep_getattr(self.encoder, name)(*x_in)
                for name, site in self.prototype_trace.iter_stochastic_nodes()
            }
        else:
            # encode with a single encoder
            res = self.encoder(*x_in)
        return res

    def forward(self, *args, **kwargs):
        """
        An automatic guide with the same ``*args, **kwargs`` as the base ``model``.

        .. note:: This method is used internally by :class:`~torch.nn.Module`.
            Users should instead use :meth:`~torch.nn.Module.__call__`.

        :return: A dict mapping sample site name to sampled value.
        :rtype: dict
        """
        # if we've never run the model before, do so now so we can inspect the model structure
        if self.prototype_trace is None:
            self._setup_prototype(*args, **kwargs)

        encoded_hidden = self.encode(*args, **kwargs)

        plates = self._create_plates(*args, **kwargs)
        result = {}
        for name, site in self.prototype_trace.iter_stochastic_nodes():
            transform = biject_to(site["fn"].support)

            with ExitStack() as stack:
                for frame in site["cond_indep_stack"]:
                    if frame.vectorized:
                        stack.enter_context(plates[frame.name])

                site_loc, site_scale = self._get_loc_and_scale(name, encoded_hidden)
                unconstrained_latent = pyro.sample(
                    name + "_unconstrained",
                    dist.Normal(
                        site_loc,
                        site_scale,
                    ).to_event(self._event_dims[name]),
                    infer={"is_auxiliary": True},
                )

                value = transform(unconstrained_latent)
                if pyro.poutine.get_mask() is False:
                    log_density = 0.0
                else:
                    log_density = transform.inv.log_abs_det_jacobian(
                        value,
                        unconstrained_latent,
                    )
                    log_density = sum_rightmost(
                        log_density,
                        log_density.dim() - value.dim() + site["fn"].event_dim,
                    )
                delta_dist = dist.Delta(
                    value,
                    log_density=log_density,
                    event_dim=site["fn"].event_dim,
                )

                result[name] = pyro.sample(name, delta_dist)

        return result

    @torch.no_grad()
    def median(self, *args, **kwargs):
        """
        Returns the posterior median value of each latent variable.

        :return: A dict mapping sample site name to median tensor.
        :rtype: dict
        """

        encoded_latent = self.encode(*args, **kwargs)

        medians = {}
        for name, site in self.prototype_trace.iter_stochastic_nodes():
            site_loc, _ = self._get_loc_and_scale(name, encoded_latent)
            median = biject_to(site["fn"].support)(site_loc)
            if median is site_loc:
                median = median.clone()
            medians[name] = median

        return medians

    @torch.no_grad()
    def quantiles(self, quantiles, *args, **kwargs):
        """
        Returns posterior quantiles each latent variable. Example::

            print(guide.quantiles([0.05, 0.5, 0.95]))

        :param quantiles: A list of requested quantiles between 0 and 1.
        :type quantiles: torch.Tensor or list
        :return: A dict mapping sample site name to a list of quantile values.
        :rtype: dict
        """

        encoded_latent = self.encode(*args, **kwargs)

        results = {}

        for name, site in self.prototype_trace.iter_stochastic_nodes():
            site_loc, site_scale = self._get_loc_and_scale(name, encoded_latent)

            site_quantiles = torch.tensor(
                quantiles, dtype=site_loc.dtype, device=site_loc.device
            )
            site_quantiles_values = dist.Normal(site_loc, site_scale).icdf(
                site_quantiles
            )
            constrained_site_quantiles = biject_to(site["fn"].support)(
                site_quantiles_values
            )
            results[name] = constrained_site_quantiles

        return results
