from contextlib import ExitStack  # python 3
from copy import deepcopy

import numpy as np
import pyro
import pyro.distributions as dist
import torch
from pyro.distributions.transforms import SoftplusTransform
from pyro.distributions.util import sum_rightmost
from pyro.infer.autoguide import AutoGuide
from pyro.infer.autoguide import AutoGuideList as PyroAutoGuideList
from pyro.infer.autoguide.guides import _deep_getattr, _deep_setattr
from pyro.infer.autoguide.utils import helpful_support_errors
from pyro.nn import PyroModule, PyroParam
from pyro.nn.module import to_pyro_module_
from torch.distributions import biject_to

from scvi._compat import Literal
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


class AutoNormalEncoder(AutoGuide):
    """
    AutoNormal posterior approximation for amortised inference,
      where mean and sd of the posterior distributions are approximated using a neural network:
      mean, sd = encoderNN(input data).

    The class supports single encoder for all parameters as well as one encoder per parameter.
    The output of encoder network is treated as a hidden layer, mean and sd are a linear function of hidden layer nodes,
    sd is transformed to positive scale using softplus. Data is log-transformed on input.

    This class requires `amortised_plate_sites` dictionary with details about amortised variables (see below).

    Guide will have the same call signature as the model, so any argument to the model can be used for encoding as
    annotated in `amortised_plate_sites`, but it does not have to be the same as observed data in the model.
    """

    def __init__(
        self,
        model,
        amortised_plate_sites: dict,
        n_in: int,
        n_hidden: dict = None,
        init_param=0,
        init_param_scale: float = 1 / 50,
        scales_offset: float = -2,
        encoder_class=FCLayersPyro,
        encoder_kwargs=None,
        multi_encoder_kwargs=None,
        encoder_instance: torch.nn.Module = None,
        create_plates=None,
        encoder_mode: Literal["single", "multiple", "single-multiple"] = "single",
    ):
        """

        Parameters
        ----------
        model
            Pyro model
        amortised_plate_sites
            Dictionary with amortised plate details:
             the name of observation/minibatch plate,
             indexes of model args to provide to encoder,
             variable names that belong to the observation plate
             and the number of dimensions in non-plate axis of each variable - such as:
             {
                 "name": "obs_plate",
                 "input": [0],  # expression data + (optional) batch index ([0, 2])
                 "input_transform": [torch.log1p], # how to transform input data before passing to NN
                 "sites": {
                     "n_s_cells_per_location": 1,
                     "y_s_groups_per_location": 1,
                     "z_sr_groups_factors": self.n_groups,
                     "w_sf": self.n_factors,
                     "l_s_add": 1,
                 }
             }
        n_in
            Number of input dimensions (for encoder_class).
        n_hidden
            Number of hidden nodes in each layer, including final layer.
        init_param
            Not implemented yet - initial values for amortised variables.
        init_param_scale
            How to scale/normalise initial values for weights converting hidden layers to mean and sd.
        encoder_class
            Class for defining encoder network.
        encoder_kwargs
            Keyword arguments for encoder_class.
        multi_encoder_kwargs
            Optional separate keyword arguments for encoder_class, useful when encoder_mode == "single-multiple".
        encoder_instance
            Encoder network instance, overrides class input and the input instance is copied with deepcopy.
        create_plates
            Function for creating plates
        encoder_mode
            Use single encoder for all variables ("single"), one encoder per variable ("multiple")
            or a single encoder in the first step and multiple encoders in the second step ("single-multiple").
        """

        super().__init__(model, create_plates=create_plates)
        self.amortised_plate_sites = amortised_plate_sites
        self.encoder_mode = encoder_mode
        self.scales_offset = scales_offset

        self.softplus = SoftplusTransform()

        if n_hidden is None:
            n_hidden = {"single": 200, "multiple": 200}
        else:
            if isinstance(n_hidden, int):
                n_hidden = {"single": n_hidden, "multiple": n_hidden}
            elif not isinstance(n_hidden, dict):
                raise ValueError("n_hidden must be either in or dict")

        encoder_kwargs = encoder_kwargs if isinstance(encoder_kwargs, dict) else dict()
        encoder_kwargs["n_hidden"] = n_hidden["single"]
        self.encoder_kwargs = encoder_kwargs
        if multi_encoder_kwargs is None:
            multi_encoder_kwargs = deepcopy(encoder_kwargs)
        self.multi_encoder_kwargs = multi_encoder_kwargs
        if "multiple" in n_hidden.keys():
            self.multi_encoder_kwargs["n_hidden"] = n_hidden["multiple"]

        self.single_n_in = n_in
        self.multiple_n_in = n_in
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
        self.encoder_instance = encoder_instance
        if "single" in self.encoder_mode:
            # create a single encoder NN
            if encoder_instance is not None:
                self.one_encoder = deepcopy(encoder_instance)
                # convert to pyro module
                to_pyro_module_(self.one_encoder)
            else:
                self.one_encoder = encoder_class(
                    n_in=self.single_n_in,
                    n_out=self.n_hidden["single"],
                    **self.encoder_kwargs
                )
            if "multiple" in self.encoder_mode:
                self.multiple_n_in = self.n_hidden["single"]

        self.init_param_scale = init_param_scale

    def _setup_prototype(self, *args, **kwargs):

        super()._setup_prototype(*args, **kwargs)

        self._event_dims = {}
        self._cond_indep_stacks = {}
        self.hidden2locs = PyroModule()
        self.hidden2scales = PyroModule()

        if "multiple" in self.encoder_mode:
            # create module for collecting multiple encoder NN
            self.multiple_encoders = PyroModule()

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

            # determine the number of hidden layers
            if "multiple" in self.encoder_mode:
                if "multiple" in self.n_hidden.keys():
                    n_hidden = self.n_hidden["multiple"]
                else:
                    n_hidden = self.n_hidden[name]
            elif "single" in self.encoder_mode:
                n_hidden = self.n_hidden["single"]
            # add linear layer for locs and scales
            param_dim = (n_hidden, self.amortised_plate_sites["sites"][name])
            init_param = np.random.normal(
                np.zeros(param_dim),
                (np.ones(param_dim) * self.init_param_scale) / np.sqrt(n_hidden),
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
                (np.ones(param_dim) * self.init_param_scale) / np.sqrt(n_hidden),
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

            if "multiple" in self.encoder_mode:
                # create multiple encoders
                if self.encoder_instance is not None:
                    # copy instances
                    encoder_ = deepcopy(self.encoder_instance).to(site["value"].device)
                    # convert to pyro module
                    to_pyro_module_(encoder_)
                    _deep_setattr(
                        self.multiple_encoders,
                        name,
                        encoder_,
                    )
                else:
                    # create instances
                    _deep_setattr(
                        self.multiple_encoders,
                        name,
                        self.encoder_class(
                            n_in=self.multiple_n_in,
                            n_out=n_hidden,
                            **self.multi_encoder_kwargs
                        ).to(site["value"].device),
                    )

    def _get_loc_and_scale(self, name, encoded_hidden):
        """
        Get mean (loc) and sd (scale) of the posterior distribution, as a linear function of encoder hidden layer.
        Parameters
        ----------
        name
            variable name
        encoded_hidden
            tensor when `encoder_mode == "single"`
            and dictionary of tensors for each site when `encoder_mode == "multiple"`

        """

        linear_locs = _deep_getattr(self.hidden2locs, name)
        linear_scales = _deep_getattr(self.hidden2scales, name)

        if "multiple" in self.encoder_mode:
            # when using multiple encoders extract hidden layer for this parameter
            encoded_hidden = encoded_hidden[name]

        locs = encoded_hidden @ linear_locs
        scales = self.softplus((encoded_hidden @ linear_scales) - self.scales_offset)

        return locs, scales

    def encode(self, *args, **kwargs):
        """
        Apply encoder network to input data to obtain hidden layer encoding.
        Parameters
        ----------
        args
            Pyro model args
        kwargs
            Pyro model kwargs
        -------

        """
        in_names = self.amortised_plate_sites["input"]
        x_in = [kwargs[i] if i in kwargs.keys() else args[i] for i in in_names]
        # apply data transform before passing to NN
        in_transforms = self.amortised_plate_sites["input_transform"]
        x_in = [in_transforms[i](x) for i, x in enumerate(x_in)]
        if "single" in self.encoder_mode:
            # encode with a single encoder
            res = self.one_encoder(*x_in)
            if "multiple" in self.encoder_mode:
                # when there is a second layer of multiple encoders fetch encoders and encode data
                x_in[0] = res
                res = {
                    name: _deep_getattr(self.multiple_encoders, name)(*x_in)
                    for name, site in self.prototype_trace.iter_stochastic_nodes()
                }
        else:
            # when there are multiple encoders fetch encoders and encode data
            res = {
                name: _deep_getattr(self.multiple_encoders, name)(*x_in)
                for name, site in self.prototype_trace.iter_stochastic_nodes()
            }
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
