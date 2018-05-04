import os
import pickle
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


class pickle_result:
    def __init__(self, filename):
        self.filename = filename

    def __call__(self, function):
        def wrapper(*args, **kwargs):
            if os.path.exists(self.filename):
                return pickle.load(open(self.filename, 'rb'))
            else:
                result = function(*args, **kwargs)
                pickle.dump(result, open(self.filename, 'wb'))
            return result

        return wrapper


def to_cuda(tensors, async=True):
    return [t.cuda(async=async) for t in tensors]
