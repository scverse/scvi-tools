import torch
from torch.autograd import Variable

from scvi.log_likelihood import compute_log_likelihood


def train(vae, data_loader_train, data_loader_test, n_epochs=20, learning_rate=0.001, kl=None):
    # Defining the optimizer
    optimizer = torch.optim.Adam(vae.parameters(), lr=learning_rate, eps=0.01)

    # Training the model
    for epoch in range(n_epochs):
        for i_batch, (sample_batch, local_l_mean, local_l_var, batch_index) in enumerate(data_loader_train):

            sample_batch = Variable(sample_batch)
            local_l_mean = Variable(local_l_mean)
            local_l_var = Variable(local_l_var)

            if torch.cuda.is_available():
                sample_batch = sample_batch.cuda(async=True)
                local_l_mean = local_l_mean.cuda(async=True)
                local_l_var = local_l_var.cuda(async=True)
                batch_index = batch_index.cuda(async=True)

            if kl is None:
                kl_ponderation = min(1, epoch / 400.)
            else:
                kl_ponderation = kl

            # Train loss is actually different from the real loss due to kl_ponderation
            if vae.batch:
                train_loss, reconst_loss, kl_divergence = vae.loss(
                    sample_batch, local_l_mean, local_l_var, kl_ponderation, batch_index)
            else:
                train_loss, reconst_loss, kl_divergence = vae.loss(
                    sample_batch, local_l_mean, local_l_var, kl_ponderation)
            real_loss = reconst_loss + kl_divergence
            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

            # Simply printing the results
        if epoch % 10 == 0:
            vae.eval()
            log_likelihood_train = compute_log_likelihood(vae, data_loader_train)
            log_likelihood_test = compute_log_likelihood(vae, data_loader_test)
            print("Epoch[%d/%d], LL-Train %.4f, LL-Test %.4f, Total Loss: %.4f, "
                  % (epoch + 1, n_epochs, log_likelihood_train, log_likelihood_test, real_loss.data[0]))
            vae.train()
