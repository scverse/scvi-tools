from contextlib import ContextDecorator

import torch
import torch.nn as nn


class enable_grad(torch.enable_grad, ContextDecorator):
    pass


class no_grad(torch.no_grad, ContextDecorator):
    pass


class eval_modules:
    def __call__(self, function):
        def wrapper(*args, **kwargs):
            [a.eval() for a in args + tuple(kwargs.items()) if isinstance(a, nn.Module)]
            result = function(*args, **kwargs)
            [a.train() for a in args + tuple(kwargs.items()) if isinstance(a, nn.Module)]
            return result

        return wrapper


def to_cuda(tensors, use_cuda=True, async=True):
    return tuple(t.cuda(async=use_cuda and async) if use_cuda else t for t in tensors)
