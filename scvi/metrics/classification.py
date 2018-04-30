import torch

from scvi.utils import no_grad, eval_modules, to_cuda


@no_grad()
@eval_modules()
def compute_accuracy(vae, data_loader, classifier=None):
    all_y_pred = []
    all_labels = []

    for i_batch, tensors in enumerate(data_loader):
        if vae.use_cuda:
            tensors = to_cuda(tensors)
        sample_batch, _, _, _, labels = tensors
        sample_batch = sample_batch.type(torch.float32)
        all_labels += [labels.view(-1)]

        if classifier is not None:
            # Then we use the specified classifier
            mu_z, _, _ = vae.z_encoder(sample_batch)
            y_pred = classifier(mu_z).argmax(dim=-1)
        else:
            # Then the vae must implement a classify function
            y_pred = vae.classify(sample_batch).argmax(dim=-1)
        all_y_pred += [y_pred]

    accuracy = (torch.cat(all_y_pred) == torch.cat(all_labels)).type(torch.float32).mean().item()

    return accuracy
