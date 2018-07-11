import torch
from scvi.utils import enable_grad


@enable_grad()
def adapt_encoder(infer, n_path=10, n_epochs=50, frequency=5):
    vae = infer.model
    parameters = list(vae.z_encoder.parameters()) + list(vae.l_encoder.parameters())
    z_encoder_state = vae.z_encoder.state_dict()
    l_encoder_state = vae.l_encoder.state_dict()
    infer.optimizer = torch.optim.Adam(parameters)
    infer.data_loaders.data_loaders_loop = [infer.data_loaders['test']]
    infer.data_loaders.to_monitor = ['test']
    infer.frequency = frequency

    # Training the model
    for i in range(n_path):
        # Re-initialize to create new path
        vae.z_encoder.load_state_dict(z_encoder_state)
        vae.l_encoder.load_state_dict(l_encoder_state)
        infer.fit(n_epochs)
    return min(infer.history["ll_test"])
