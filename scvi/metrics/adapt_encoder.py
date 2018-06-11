import torch

from scvi.metrics.stats import Stats
from scvi.utils import to_cuda, enable_grad


@enable_grad()
def adapt_encoder(vae, dataloader, n_path=10, n_epochs=50, record_freq=5):
    parameters = list(vae.z_encoder.parameters()) + list(vae.l_encoder.parameters())
    z_encoder_state = vae.z_encoder.state_dict()
    l_encoder_state = vae.l_encoder.state_dict()

    optimizer = torch.optim.Adam(parameters, eps=0.01)
    # Getting access to the stats during training
    stats = Stats(n_epochs=n_epochs, record_freq=record_freq, names=['test'], verbose=False)
    stats.callback(vae, dataloader)
    best_ll = stats.history["LL_test"][0]

    # Training the model
    for i in range(n_path):
        # Re-initialize to create new path
        vae.z_encoder.load_state_dict(z_encoder_state)
        vae.l_encoder.load_state_dict(l_encoder_state)
        for epoch in range(n_epochs):
            for i_batch, tensors in enumerate(dataloader):
                if vae.use_cuda:
                    tensors = to_cuda(tensors)
                sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
                sample_batch = sample_batch.type(torch.float32)

                reconst_loss, _ = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index, y=labels)
                train_loss = torch.mean(reconst_loss)

                optimizer.zero_grad()
                train_loss.backward()
                optimizer.step()

            stats.callback(vae, dataloader)
        best_ll = min(min(stats.history["LL_test"]), best_ll)
    return best_ll
