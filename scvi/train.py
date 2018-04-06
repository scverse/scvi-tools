import numpy as np
import torch
from torch.autograd import Variable
from torch.utils.data import DataLoader

from scvi.log_likelihood import log_zinb_positive

if torch.cuda.is_available():
    dtype = torch.cuda.FloatTensor
else:
    dtype = torch.FloatTensor


def train(vae, data_loader, n_epochs=20, learning_rate=0.001, kl=None):
    # Defining the optimizer
    optimizer = torch.optim.Adam(vae.parameters(), lr=learning_rate, eps=0.01)
    iter_per_epoch = len(data_loader)

    # Training the model
    for epoch in range(n_epochs):
        for i_batch, (sample_batched, local_l_mean, local_l_var, batch_index) in enumerate(data_loader):
            sample_batched = Variable(sample_batched, requires_grad=False)
            local_l_mean = Variable(local_l_mean.type(dtype), requires_grad=False)
            local_l_var = Variable(local_l_var.type(dtype), requires_grad=False)

            px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v = vae(sample_batched)

            # Computing the reconstruction loss

            reconst_loss = -log_zinb_positive(sample_batched, px_rate, torch.exp(px_r), px_dropout)

            # Computing the kl divergence
            kl_divergence_z = torch.sum(0.5 * (qz_m ** 2 + qz_v - torch.log(qz_v + 1e-8) - 1), dim=1)
            kl_divergence_l = torch.sum(0.5 * (
                ((ql_m - local_l_mean) ** 2) / local_l_var + ql_v / local_l_var
                + torch.log(local_l_var + 1e-8) - torch.log(ql_v + 1e-8) - 1), dim=1)

            if kl is None:
                kl_ponderation = np.minimum(1, epoch / 400.)
            else:
                kl_ponderation = kl

            kl_ponderation = Variable(torch.from_numpy(np.array([kl_ponderation])).type(dtype),
                                      requires_grad=False)
            kl_divergence = (kl_divergence_z + kl_divergence_l)

            # Backprop + Optimize
            total_loss = torch.mean(reconst_loss + kl_ponderation * kl_divergence)
            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            # Simply printing the results
            if i_batch % 100 == 0 and epoch % 1 == 0:
                print("Epoch[%d/%d], Step [%d/%d], Total Loss: %.4f, "
                      "Reconst Loss: %.4f, KL Div: %.7f"
                      % (epoch + 1, n_epochs, i_batch + 1, iter_per_epoch, total_loss.data[0],
                         torch.mean(reconst_loss).data[0], torch.mean(kl_divergence).data[0]))


def compute_log_likelihood(vae, gene_dataset):
    data_loader_test = DataLoader(gene_dataset, batch_size=gene_dataset.total_size, shuffle=False, num_workers=1)
    for i_batch, (sample_batched, local_l_mean, local_l_var, batch_index) in enumerate(data_loader_test):
        sample_batched = Variable(sample_batched, requires_grad=False)
        px_scale, px_r, px_rate, px_dropout, qz_m, qz_v, ql_m, ql_v = vae(sample_batched)
        log_lkl = torch.mean(-log_zinb_positive(sample_batched, px_rate, torch.exp(px_r), px_dropout)).data[
            0]
    return log_lkl
