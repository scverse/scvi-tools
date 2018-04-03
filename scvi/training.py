import torch
import numpy as np
from torch.utils.data import DataLoader
import scvi.dataset as dt
import scvi.log_likelihood as lkl
from torch.autograd import Variable

# Generating samples according to a ZINB process
batch_size = 20
nb_genes = 100
data = np.random.negative_binomial(5, 0.3, size=(batch_size, nb_genes))
newdata = np.ones((batch_size, nb_genes))
mask = np.random.binomial(n=1, p=0.7, size=(batch_size, nb_genes))
for i in range(batch_size):
    newdata[i, :] = data[i, :] / np.sum(data[i, :])
    newdata[i, :] = newdata[i, :] * mask[i, :]
data = torch.Tensor(newdata)


def train(vae, data=data, n_epochs=20):
    # Defining the optimizer
    optimizer = torch.optim.Adam(vae.parameters(), lr=0.001)

    # For now, just numeric mode. The values chosen are just exemples.
    library_size_mean, library_size_var = np.array([10]), np.array([3])
    local_l_mean = Variable((torch.from_numpy(library_size_mean)).type(torch.FloatTensor))
    local_l_var = Variable((torch.from_numpy(library_size_var)).type(torch.FloatTensor))

    # Creating a GeneExpressionDataset and a DataLoader
    gene_dataset = dt.GeneExpressionDataset(data)
    data_loader = DataLoader(gene_dataset, batch_size=4,
                             shuffle=True, num_workers=4)

    iter_per_epoch = len(data_loader)

    # Training the model
    for epoch in range(n_epochs):
        for i_batch, sample_batched in enumerate(data_loader):
            sample_batched = Variable(sample_batched)

            out, px_r, px_rate, px_dropout, mu_z, log_var_z, mu_l, log_var_l = vae(sample_batched)

            # Checking if the approximation for the log-likelihood is valid
            if np.abs(lkl.log_zinb_positive_approx(sample_batched, px_rate, torch.exp(px_r), px_dropout).data[0] -
                      lkl.log_zinb_positive_real(sample_batched, px_rate, torch.exp(px_r), px_dropout).data[0]) > 1e-1:
                print("Error due to numerical approximation is: ")
                print(
                    np.abs(lkl.log_zinb_positive_approx(sample_batched, px_rate, torch.exp(px_r), px_dropout).data[0] -
                           lkl.log_zinb_positive_real(sample_batched, px_rate, torch.exp(px_r), px_dropout).data[0]))
                print("Careful with the numerical approximations!")

            # Computing the reconstruction loss
            reconst_loss = -lkl.log_zinb_positive_approx(sample_batched, px_rate, torch.exp(px_r), px_dropout)

            # Computing the kl divergence
            kl_divergence_z = torch.sum(0.5 * (mu_z ** 2 + torch.exp(log_var_z) - log_var_z - 1))
            kl_divergence_l = torch.sum(0.5 * (
                ((mu_l - local_l_mean) ** 2) / local_l_var + torch.exp(log_var_l) / local_l_var
                + torch.log(local_l_var) - log_var_l - 1))
            kl_divergence = kl_divergence_z + kl_divergence_l

            # Backprop + Optimize
            total_loss = reconst_loss + kl_divergence
            optimizer.zero_grad()
            total_loss.backward()
            optimizer.step()

            # Simply printing the results
            if i_batch % 100 == 0 and epoch % 1 == 0:
                print("Epoch[%d/%d], Step [%d/%d], Total Loss: %.4f, "
                      "Reconst Loss: %.4f, KL Div: %.7f"
                      % (epoch + 1, 20, i_batch + 1, iter_per_epoch, total_loss.data[0],
                         reconst_loss.data[0], kl_divergence.data[0]))
