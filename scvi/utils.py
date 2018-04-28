import torch


def one_hot(index, n_cat):
    onehot = torch.zeros(index.size(0), n_cat, device=index.device)
    onehot.scatter_(1, index, 1)
    return onehot.type(torch.float32)


def enumerate_discrete(x, y_dim):
    def batch(batch_size, label):
        labels = (torch.ones(batch_size, 1, device=x.device, dtype=torch.long) * label)
        return one_hot(labels, y_dim)

    batch_size = x.size(0)
    return torch.cat([batch(batch_size, i) for i in range(y_dim)])


def to_cuda(tensor_list, async=True):
    return [t.cuda(async=async) for t in tensor_list]


def compute_accuracy(vae, data_loader, classifier=None):
    all_y_pred = []
    all_labels = []

    with torch.no_grad():
        for i_batch, (sample_batch, _, _, _, labels) in enumerate(data_loader):
            sample_batch = sample_batch.type(torch.float32)
            if vae.using_cuda:
                sample_batch = sample_batch.cuda(async=True)
                labels = labels.cuda(async=True)
            all_labels += [labels.view(-1)]

            if classifier is not None:
                # Then we use the specified classifier
                mu_z, _, _ = vae.z_encoder(sample_batch)
                y_pred = classifier(mu_z).argmax(dim=-1)
            else:
                # Then the vae must implement a classify function
                y_pred = vae.classify(sample_batch).argmax(dim=-1)
            all_y_pred += [y_pred]

    accuracy = (torch.cat(all_y_pred) == torch.cat(all_labels)).type(torch.FloatTensor).mean().item()

    return accuracy
