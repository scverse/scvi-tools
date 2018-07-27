import copy


def adapt_encoder(infer, n_path=10, n_epochs=50, frequency=5):
    vae = infer.model
    params = list(vae.z_encoder.parameters()) + list(vae.l_encoder.parameters())
    z_encoder_state = copy.deepcopy(vae.z_encoder.state_dict())
    l_encoder_state = copy.deepcopy(vae.l_encoder.state_dict())
    infer.data_loaders.loop = ['test']
    infer.data_loaders.to_monitor = ['test']
    infer.frequency = frequency

    # Training the model
    for i in range(n_path):
        # Re-initialize to create new path
        vae.z_encoder.load_state_dict(z_encoder_state)
        vae.l_encoder.load_state_dict(l_encoder_state)
        infer.train(n_epochs, params=params)
    return min(infer.history["ll_test"])
