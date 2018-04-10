import numpy as np
import torch
from torch.autograd import Variable

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
        for i_batch, (sample_batch, local_l_mean, local_l_var, batch_index) in enumerate(data_loader):
            sample_batch = Variable(sample_batch.type(dtype), requires_grad=False)
            local_l_mean = Variable(local_l_mean.type(dtype), requires_grad=False)
            local_l_var = Variable(local_l_var.type(dtype), requires_grad=False)

            if kl is None:
                kl_ponderation = np.minimum(1, epoch / 400.)
            else:
                kl_ponderation = kl

            # Train loss is actually different from the real loss due to kl_ponderation
            train_loss, reconst_loss, kl_divergence = vae.loss(sample_batch, local_l_mean, local_l_var, kl_ponderation)
            real_loss = reconst_loss + kl_divergence
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

            # Simply printing the results
            if i_batch % 100 == 0 and epoch % 1 == 0:
                print("Epoch[%d/%d], Step [%d/%d], Total Loss: %.4f, "
                      "Reconst Loss: %.4f, KL Div: %.7f"
                      % (epoch + 1, n_epochs, i_batch + 1, iter_per_epoch, real_loss.data[0],
                         torch.mean(reconst_loss).data[0], torch.mean(kl_divergence).data[0]))
