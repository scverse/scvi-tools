import torch
from torch.distributions.constraints import Constraint


class _Optional(Constraint):
    def __init__(self, constraint: Constraint):
        self.constraint = constraint

    def check(self, value: torch.Tensor) -> torch.Tensor:
        if value is None:
            return torch.ones(1, dtype=torch.bool)
        return self.constraint.check(value)

    def __repr__(self) -> str:
        return f"Optional({self.constraint})"


def optional_constraint(constraint: Constraint) -> Constraint:
    """Returns a wrapped constraint that allows optional values."""
    return _Optional(constraint)
