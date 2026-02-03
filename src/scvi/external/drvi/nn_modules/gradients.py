import torch


class GradScale(torch.autograd.Function):
    @staticmethod
    def forward(ctx, x, scale):
        ctx.scale = scale
        return x  # forward pass unchanged

    @staticmethod
    def backward(ctx, grad_output):
        return grad_output * ctx.scale, None  # scale gradient only


def grad_scale(x, scale):
    return GradScale.apply(x, scale)


class GradientScaler(torch.nn.Module):
    def __init__(self, scale: float):
        super().__init__()
        self.register_buffer("scale", torch.tensor(scale, dtype=torch.float32))

    def forward(self, x):
        return grad_scale(x, self.scale)

    def extra_repr(self):
        return f"scale={self.scale.item()}"
