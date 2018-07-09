import torch

from scvi.metrics.early_stopping import InferenceTask
from scvi.utils import to_cuda, enable_grad


@enable_grad()
def adapt_encoder(infer, data_loader, n_path=10, n_epochs=50, record_freq=5, use_cuda=True):
    vae = infer.model
    parameters = list(vae.z_encoder.parameters()) + list(vae.l_encoder.parameters())
    z_encoder_state = vae.z_encoder.state_dict()
    l_encoder_state = vae.l_encoder.state_dict()

    optimizer = torch.optim.Adam(parameters, eps=0.01)
    # Getting access to the stats during training
    stats = InferenceTask(n_epochs=n_epochs, record_freq=record_freq, names=['test'], verbose=False, use_cuda=use_cuda)
    stats.callback(vae, data_loader)
    best_ll = stats.history["LL_test"][0]

    # Training the model
    for i in range(n_path):
        # Re-initialize to create new path
        vae.z_encoder.load_state_dict(z_encoder_state)
        vae.l_encoder.load_state_dict(l_encoder_state)
        for epoch in range(n_epochs):
            for i_batch, tensors in enumerate(data_loader):
                tensors = to_cuda(tensors, use_cuda=use_cuda)
                sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors
                sample_batch = sample_batch.type(torch.float32)

                reconst_loss, _ = vae(sample_batch, local_l_mean, local_l_var, batch_index=batch_index, y=labels)
                train_loss = torch.mean(reconst_loss)

                optimizer.zero_grad()
                train_loss.backward()
                optimizer.step()

            stats.callback(vae, data_loader)
        best_ll = min(min(stats.history["LL_test"]), best_ll)
    return best_ll
