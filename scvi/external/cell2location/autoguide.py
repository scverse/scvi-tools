from contextlib import ExitStack  # python 3

import numpy as np
import pyro
import pyro.distributions as dist
import torch
from pyro.distributions import constraints
from pyro.distributions.transforms import SoftplusTransform
from pyro.distributions.util import sum_rightmost
from pyro.infer.autoguide import AutoGuide
from pyro.infer.autoguide.utils import helpful_support_errors
from pyro.nn import PyroModule
from torch.distributions import biject_to, transform_to

from scvi.nn import FCLayers


class FCLayersPyro(FCLayers, PyroModule):
    pass


@biject_to.register(constraints.positive)
@transform_to.register(constraints.positive)
def _transform_to_positive(constraint):
    return SoftplusTransform()


# sample global parameters automatically
#        auto_guide_kwargs = auto_guide_kwargs if isinstance(auto_guide_kwargs, dict) else dict()
#        self.auto_guide = AutoNormal(
#            pyro.poutine.block(model, hide=list(amortised_plate_sites['sites'].keys())),
#            **auto_guide_kwargs
#        )


class AutoNormalEncoder(AutoGuide):
    def __init__(
        self,
        model,
        amortised_plate_sites,
        n_in,
        encoder_kwargs=None,
        create_plates=None,
    ):

        encoder_kwargs = encoder_kwargs if isinstance(encoder_kwargs, dict) else dict()

        super().__init__(model, create_plates=create_plates)
        self.amortised_plate_sites = amortised_plate_sites

        self.softplus = SoftplusTransform()
        # create encoder NN
        n_out = (
            np.sum(
                [
                    np.sum(amortised_plate_sites["sites"][k])
                    for k in amortised_plate_sites["sites"].keys()
                ]
            )
            * 2
        )
        self.encoder = FCLayersPyro(n_in=n_in, n_out=n_out, **encoder_kwargs)
        # create indices for loc and scales of each site
        counter = 0
        self.indices = dict()
        for site, n_dim in amortised_plate_sites["sites"].items():
            self.indices[site] = {
                "locs": np.arange(counter, counter + n_dim),
                "scales": np.arange(counter + n_dim, counter + n_dim * 2),
            }
            counter += n_dim * 2

    def _setup_prototype(self, *args, **kwargs):

        super()._setup_prototype(*args, **kwargs)

        self._event_dims = {}
        self._cond_indep_stacks = {}

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

    def _get_loc_and_scale(self, name, encoded_latent):

        ind_locs = (
            torch.tensor(self.indices[name]["locs"], device=encoded_latent.device)
            .long()
            .squeeze()
        )
        ind_scales = (
            torch.tensor(self.indices[name]["scales"], device=encoded_latent.device)
            .long()
            .squeeze()
        )

        locs = torch.index_select(encoded_latent, -1, ind_locs)
        scales = self.softplus(torch.index_select(encoded_latent, -1, ind_scales) - 2)

        return locs, scales

    def encode(self, *args, **kwargs):
        in_names = self.amortised_plate_sites["in"]
        x_in = [
            kwargs[in_name] if in_name in kwargs.keys() else args[i]
            for i, in_name in enumerate(in_names)
        ]
        # apply log1p
        x_in = [torch.log1p(x) for x in x_in]
        return self.encoder(*x_in)

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

        encoded_latent = self.encode(*args, **kwargs)

        plates = self._create_plates(*args, **kwargs)
        result = {}
        for name, site in self.prototype_trace.iter_stochastic_nodes():
            transform = biject_to(site["fn"].support)

            with ExitStack() as stack:
                for frame in site["cond_indep_stack"]:
                    if frame.vectorized:
                        stack.enter_context(plates[frame.name])

                site_loc, site_scale = self._get_loc_and_scale(name, encoded_latent)
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
