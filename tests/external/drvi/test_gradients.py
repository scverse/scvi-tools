import torch

from scvi.external.drvi.nn_modules.gradients import GradientScaler


class TestGradientScaler:
    """Test cases for GradientScaler class.

    The rest of the parts are used here, so not separate test.
    """

    def test_gradient_scaler_forward_pass(self):
        """Test that forward pass returns input unchanged."""
        scaler = GradientScaler(scale=2.0)
        x = torch.randn(3, 4, requires_grad=True)

        output = scaler(x)

        # Forward pass should return input unchanged
        assert torch.allclose(output, x)
        assert output.requires_grad == x.requires_grad

    def test_gradient_scaler_gradient_scaling(self):
        """Test that gradients are properly scaled during backward pass."""
        scales_to_test = [-1.0, 0.1, 1.0, 2.5, 10.0]
        for scale_value in scales_to_test:
            scaler = GradientScaler(scale=scale_value)
            x = torch.randn(2, 3, requires_grad=True)

            # Forward pass
            output = scaler(x)
            loss = (output**2).sum()

            # Backward pass
            loss.backward()

            # Gradient should be scaled by the scale factor
            expected_grad = 2 * x * scale_value  # 2*x from loss derivative
            assert torch.allclose(x.grad, expected_grad)
