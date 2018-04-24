import random

import torch
from torch.autograd import Variable

from scvi.log_likelihood import compute_log_likelihood


def train(vae, data_loader_train, data_loader_test, n_epochs=20, learning_rate=0.001, kl=None,
          early_stopping_criterion=(20, 0.01)):
    # Defining the optimizer
    optimizer = torch.optim.Adam(vae.parameters(), lr=learning_rate, eps=0.01)

    # initialize
    (patience, threshold) = (early_stopping_criterion[0], early_stopping_criterion[1])
    current_performances_kl = torch.ones(patience)
    current_performances_reconst = torch.ones(patience)

    # Training the model
    for epoch in range(n_epochs):
        # initialize kl, reconst
        total_current_kl = 0
        total_current_reconst = 0
        for i_batch, (sample_batch, local_l_mean, local_l_var, batch_index, labels) in enumerate(data_loader_train):
            sample_batch = sample_batch.type(torch.FloatTensor)
            sample_batch = Variable(sample_batch)
            local_l_mean = Variable(local_l_mean)
            local_l_var = Variable(local_l_var)
            labels = labels if random.random() < 0.5 else None
            if vae.using_cuda:
                sample_batch = sample_batch.cuda(async=True)
                local_l_mean = local_l_mean.cuda(async=True)
                local_l_var = local_l_var.cuda(async=True)
                batch_index = batch_index.cuda(async=True)
                # labels = labels.cuda(async=True)

            if kl is None:
                kl_ponderation = min(1, epoch / 400.)
            else:
                kl_ponderation = kl

            # Train loss is actually different from the real loss due to kl_ponderation
            if vae.batch:
                reconst_loss, kl_divergence = vae(
                    sample_batch, local_l_mean, local_l_var, batch_index=batch_index, y=labels)
            else:
                reconst_loss, kl_divergence = vae(
                    sample_batch, local_l_mean, local_l_var, y=labels)
            total_current_kl += torch.mean(kl_divergence)
            total_current_reconst += torch.mean(reconst_loss)

            train_loss = torch.mean(reconst_loss) + kl_ponderation * torch.mean(kl_divergence)
            real_loss = torch.mean(reconst_loss) + torch.mean(kl_divergence)

            optimizer.zero_grad()
            train_loss.backward()
            optimizer.step()

        current_performances_kl[:-1] = current_performances_kl[1:]
        current_performances_reconst[:-1] = current_performances_reconst[1:]
        current_performances_kl[-1] = total_current_kl.data[0]
        current_performances_reconst[-1] = total_current_reconst.data[0]

        # Computing the relative improvment of kl and reconstruction loss
        # over the chosen number of epochs
        reconst_relative_improvement = ((current_performances_reconst[0] - current_performances_reconst[-1])
                                        / current_performances_reconst[0])
        kl_relative_improvement = ((current_performances_kl[0] - current_performances_kl[-1])
                                   / current_performances_kl[0])
        # Test whether stopping criterions are met
        if epoch > patience and kl_relative_improvement < threshold and reconst_relative_improvement < threshold:
            # We then stop the iterations
            print("Stopping the training after %d epochs: over the %d last epochs,"
                  " kl divergence improvement was %.4f, reconstruction loss improvement was %.4f"
                  % (epoch, patience, kl_relative_improvement, reconst_relative_improvement))

            break

            # Simply printing the results
        if epoch % 10 == 0:
            vae.eval()
            # No need to compute it twice
            log_likelihood_train = total_current_reconst.data[0] + total_current_kl.data[0]
            log_likelihood_test = compute_log_likelihood(vae, data_loader_test)
            print("Epoch[%d/%d], LL-Train %.4f, LL-Test %.4f, Total Loss: %.4f, "
                  % (epoch + 1, n_epochs, log_likelihood_train, log_likelihood_test, real_loss.data[0]))
            vae.train()
