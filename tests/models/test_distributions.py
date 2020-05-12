from unittest import TestCase
import torch
from .old_distributions import NegativeBinomial as OldNegativeBinomial
from .old_distributions import log_nb_positive as old_log_nb_positive
from scvi.models.distributions import NegativeBinomial
from scvi.models.log_likelihood import log_nb_positive


class TestNegativeBinomial(TestCase):
    def setUp(self):
        self.x = torch.Tensor([[1, 2, 3],
                               [4, 5, 6]])
        self.mu = torch.Tensor([[10, 20, 30],
                                [10, 20, 30]])
        self.theta = torch.Tensor([[2, 1, 3],
                                   [2, 1, 3]])

    def test_log_nb_positive(self):
        exp = old_log_nb_positive(self.x, self.mu, self.theta)
        res = log_nb_positive(self.x, torch.log(self.mu),
                              torch.log(self.theta))
        self.assertTrue(torch.allclose(exp, res))

    def test_log_equivalence(self):
        exp = OldNegativeBinomial(mu=self.mu, theta=self.theta)
        res = NegativeBinomial(log_mu=torch.log(self.mu), log_theta=torch.log(self.theta))
        self.assertTrue(torch.allclose(torch.log(exp.mu), res.log_mu))
        self.assertTrue(torch.allclose(torch.log(exp.theta), res.log_theta))


