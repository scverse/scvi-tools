import matplotlib.pyplot as plt
import torch

from scvi.inference import UnsupervisedTrainer

plt.switch_backend('agg')

class AutoZITrainer(UnsupervisedTrainer):

    def __init__(self, model, gene_dataset, kl_global=True, train_size=0.8, test_size=None, **kwargs):
        self.kl_global = kl_global
        super().__init__(model, gene_dataset, **kwargs)
        if type(self) is AutoZITrainer:
            self.train_set, self.test_set, self.validation_set = self.train_test_validation(
                model, gene_dataset, train_size, test_size
            )
            self.train_set.to_monitor = ["elbo"]
            self.test_set.to_monitor = ["elbo"]
            self.validation_set.to_monitor = ["elbo"]

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, y = tensors
        reconst_loss, kl_divergence, kl_divergence_global =\
            self.model(sample_batch, local_l_mean, local_l_var, batch_index, y)

        loss = len(self.train_set.indices) * torch.mean(reconst_loss + self.kl_weight * kl_divergence) \
               + self.kl_global * kl_divergence_global
        return loss
