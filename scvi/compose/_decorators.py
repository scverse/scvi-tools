from collections.abc import Mapping, Sequence
from functools import wraps
from typing import Any, Callable, Union

import torch
from torch.nn import Module


def auto_move_data(fn: Callable) -> Callable:
    """
    Decorator for :class:`~torch.nn.Module methods to move data to correct device.

    Input arguments are moved automatically to the correct device.
    It has no effect if applied to a method of an object that is not an instance of
    :class:`~torch.nn.Module` and is typically applied to ``__call__``
    or ``forward``.

    Parameters
    ----------
    fn
        A nn.Module method for which the arguments should be moved to the device
        the parameters are on.
    """

    @wraps(fn)
    def auto_transfer_args(self, *args, **kwargs):
        if not isinstance(self, Module):
            return fn(self, *args, **kwargs)

        # decorator only necessary after training
        if self.training:
            return fn(self, *args, **kwargs)

        device = list(set(p.device for p in self.parameters()))
        if len(device) > 1:
            raise RuntimeError("Model tensors on multiple devices.")
        else:
            device = device[0]
        args = _move_data_to_device(args, device)
        kwargs = _move_data_to_device(kwargs, device)
        return fn(self, *args, **kwargs)

    return auto_transfer_args


def _move_data_to_device(batch: Any, device: torch.device):
    """
    Transfers a collection of data to the given device.

    Any object that defines a method ``to(device)`` will be moved and all other objects
    in the collection will be left untouched.

    Parameters
    ----------
    batch
        A tensor or collection of tensors or anything that has a method `.to(...)`.
        See :func:`apply_to_collection` for a list of supported collection types.
    device
        The device to which the data should be moved

    Returns
    -------
        The same collection but with all contained tensors residing on the new device.
    """

    def batch_to(data):
        kwargs = dict(non_blocking=True) if isinstance(data, torch.Tensor) else {}
        return data.to(device, **kwargs)

    return _apply_to_collection(batch, dtype=torch.Tensor, function=batch_to)


def _apply_to_collection(
    data: Any, dtype: Union[type, tuple], function: Callable, *args, **kwargs
) -> Any:
    """
    Recursively applies a function to all elements of a certain dtype.

    Parameters
    ----------
    data
        The collection to apply the function to
    dtype
        The given function will be applied to all elements of this dtype
    function
        The function to apply
    *args
        positional arguments (will be forwarded to calls of ``function``)
    **kwargs
        keyword arguments (will be forwarded to calls of ``function``)

    Returns
    -------
    The resulting collection
    """
    elem_type = type(data)

    # Breaking condition
    if isinstance(data, dtype):
        return function(data, *args, **kwargs)

    # Recursively apply to collection items
    elif isinstance(data, Mapping):
        return elem_type(
            {
                k: _apply_to_collection(v, dtype, function, *args, **kwargs)
                for k, v in data.items()
            }
        )
    elif isinstance(data, tuple) and hasattr(data, "_fields"):  # named tuple
        return elem_type(
            *(_apply_to_collection(d, dtype, function, *args, **kwargs) for d in data)
        )
    elif isinstance(data, Sequence) and not isinstance(data, str):
        return elem_type(
            [_apply_to_collection(d, dtype, function, *args, **kwargs) for d in data]
        )

    # data is neither of dtype, nor a collection
    return data
