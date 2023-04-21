from typing import List

import torch


def subset_distribution(
    my_distribution: torch.distributions.Distribution,
    index: torch.Tensor,
    dim: int = 0,
) -> torch.distributions.Distribution:
    """Utility function to subset the parameter of a Pytorch distribution."""
    return my_distribution.__class__(
        **{
            name: torch.index_select(
                getattr(my_distribution, name), dim=dim, index=index
            )
            for name in my_distribution.arg_constraints.keys()
        }
    )


class DistributionsConcatenator:
    """Utility class to concatenate Pytorch distributions and move them to cpu."""

    def __init__(self):
        self.storage = {}

    def add_distributions(self, forward_outputs: dict):
        """Add a dictionary of distributions to the concatenator."""
        for key, potential_distribution in forward_outputs.items():
            if isinstance(potential_distribution, torch.distributions.Distribution):
                if key not in self.storage:
                    params = {
                        name: []
                        for name in potential_distribution.arg_constraints.keys()
                    }
                    self.storage[key] = dict(
                        cls=potential_distribution.__class__,
                        **params,
                    )
                new_params = {
                    name: getattr(potential_distribution, name).cpu()
                    for name in potential_distribution.arg_constraints.keys()
                }
                for param_name, param in new_params.items():
                    self.storage[key][param_name].append(param)

    @staticmethod
    def _find_concat_dim(my_list: List):
        ndims = my_list[0].ndim
        if ndims == 2:
            return 0
        elif ndims == 3:
            return 1
        else:
            raise ValueError("Only 2D and 3D tensors are supported.")

    def get_concatenated_distributions(self):
        """Returns concatenated distributions."""
        dists = {}
        for dist_name, dist_props in self.storage.items():
            dist_cls = dist_props.pop("cls")
            concat_params = {
                key: torch.cat(value, dim=self._find_concat_dim(value))
                for key, value in dist_props.items()
            }
            dists[dist_name] = dist_cls(**concat_params)
        return dists
