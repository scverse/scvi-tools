from __future__ import annotations

from typing import TYPE_CHECKING

from torch.distributions import Normal as NormalTorch

if TYPE_CHECKING:
    import torch


class Normal(NormalTorch):
    """Normal distribution.

    Parameters
    ----------
    mu
        loc of the Normal distribution.
    scale
        scale of the Normal distribution.
    validate_args
        whether to validate input.
    normal_mu
        Normalized mean expression of the distribution.
        This optional parameter is not used in any computations, but allows to store
        normalization expression levels.

    """

    def __init__(
        self,
        mu: torch.Tensor,
        scale: torch.Tensor,
        validate_args: bool | None = None,
        normal_mu: torch.Tensor = None,
    ):
        super().__init__(loc=mu, scale=scale, validate_args=validate_args)
        self.normal_mu = normal_mu

    def __repr__(self) -> str:
        param_names = [k for k, _ in self.arg_constraints.items() if k in self.__dict__]
        args_string = ", ".join(
            [
                f"{p}: "
                f"{self.__dict__[p] if self.__dict__[p].numel() == 1 else self.__dict__[p].size()}"
                for p in param_names
                if self.__dict__[p] is not None
            ]
        )
        return self.__class__.__name__ + "(" + args_string + ")"

    def get_normalized(self, key) -> torch.Tensor:
        if key == "mu":
            return self.loc
        elif key == "scale":
            return self.normal_mu
        else:
            raise ValueError(f"normalized key {key} not recognized")
