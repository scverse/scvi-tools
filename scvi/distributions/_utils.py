import torch


def subset_distribution(
    my_distribution: torch.distributions.Distribution,
    index: torch.Tensor,
    dim: int = 0,
) -> torch.distributions.Distribution:
    """Utility function to subset the parameter of a Pytorch distribution."""
    return my_distribution.__class__(
        **{
            name: torch.index_select(getattr(my_distribution, name), dim=dim, index=index)
            for name in my_distribution.arg_constraints.keys()
        }
    )


class DistributionConcatenator:
    """Utility class to concatenate Pytorch distributions and move them to cpu.

    All distributions must be of the same type.
    """

    def __init__(self):
        self._params = None
        self.distribution_cls = None

    def store_distribution(self, dist: torch.distributions.Distribution):
        """Add a dictionary of distributions to the concatenator.

        Parameters
        ----------
        dist:
            A Pytorch distribution.
        """
        if self._params is None:
            self._params = {name: [] for name in dist.arg_constraints.keys()}
            self.distribution_cls = dist.__class__
        new_params = {name: getattr(dist, name).cpu() for name in dist.arg_constraints.keys()}
        for param_name, param in new_params.items():
            self._params[param_name].append(param)

    def get_concatenated_distributions(self, axis=0):
        """Returns a concatenated `Distribution` object along the specified axis."""
        concat_params = {key: torch.cat(value, dim=axis) for key, value in self._params.items()}
        return self.distribution_cls(**concat_params)
