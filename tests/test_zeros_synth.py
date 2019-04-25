import numpy as np
import torch

from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
from scvi.dataset.synthetic import ZISyntheticDatasetCorr
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def test_zeros_classif():
    synth_data = ZISyntheticDatasetCorr(lam_dropout=0.1, dropout_coef=0.5)
    zeros_mask = (synth_data.X == 0).astype(bool)
    is_technical_all = synth_data.is_technical[0, :]
    is_technical_gt = is_technical_all[zeros_mask]  # 1d array

    mdl = VAE(n_input=synth_data.nb_genes, n_batch=synth_data.n_batches,
              reconstruction_loss='zinb')
    trainer = UnsupervisedTrainer(model=mdl, gene_dataset=synth_data, use_cuda=True, train_size=0.8,
                                  frequency=1,
                                  early_stopping_kwargs={
                                     'early_stopping_metric': 'll',
                                     'save_best_state_metric': 'll',
                                     'patience': 15,
                                     'threshold': 3,
    })
    trainer.train(n_epochs=150, lr=1e-3)
    full = trainer.create_posterior(trainer.model, synth_data,
                                    indices=np.arange(len(synth_data)))
    # Infer if cell/gene zeros have technical origin
    is_technical_infer = []

    with torch.no_grad():
        for tensors in full.sequential():
            sample_batch, _, _, batch_index, labels = tensors
            px_scale, px_dispersion, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v, library = mdl.inference(
                sample_batch, batch_index)

            is_technical_batch = torch.zeros((sample_batch.size(0), sample_batch.size(1), 100))
            for n_mc_sim in range(100):
                p = px_rate / (px_rate + px_dispersion)
                r = px_dispersion
                l_train = torch.distributions.Gamma(concentration=r, rate=(1 - p) / p).sample()
                l_train = torch.clamp(l_train, max=1e18)
                X = torch.distributions.Poisson(l_train).sample()
                p_zero = 1.0 / (1.0 + torch.exp(-px_dropout))
                random_prob = torch.rand_like(p_zero)
                X[random_prob <= p_zero] = 0
                is_technical_batch[:, :, n_mc_sim] = (random_prob <= p_zero)
                print(torch.mean(is_technical_batch, dim=(-1)))
            is_technical_batch = torch.mean(is_technical_batch, dim=(-1)) >= 0.5
            print(px_dropout.min(), px_dropout.max())
            is_technical_infer.append(is_technical_batch.cpu().numpy())
    is_technical_infer_all = np.concatenate(is_technical_infer)
    is_technical_infer_zeros = is_technical_infer_all[zeros_mask]
    assert is_technical_infer_all.shape == zeros_mask.shape
    assert is_technical_infer_zeros.shape == is_technical_gt.shape

    # confusion_mat = confusion_matrix(is_technical_gt, is_technical_infer_zeros, labels=[0, 1])

    def plot_confusion_matrix(y_true, y_pred, classes=[0, 1],
                              normalize=False,
                              title=None,
                              cmap=plt.cm.Blues):
        """
        This function prints and plots the confusion matrix.
        Normalization can be applied by setting `normalize=True`.
        """
        if not title:
            if normalize:
                title = 'Normalized confusion matrix'
            else:
                title = 'Confusion matrix, without normalization'

        # Compute confusion matrix
        cm = confusion_matrix(y_true, y_pred)
        # Only use the labels that appear in the data
        if normalize:
            cm = cm.astype('float') / cm.sum(axis=1)[:, np.newaxis]
            print("Normalized confusion matrix")
        else:
            print('Confusion matrix, without normalization')
        print(cm)

        fig, ax = plt.subplots()
        im = ax.imshow(cm, interpolation='nearest', cmap=cmap)
        ax.figure.colorbar(im, ax=ax)
        # We want to show all ticks...
        ax.set(xticks=np.arange(cm.shape[1]),
               yticks=np.arange(cm.shape[0]),
               # ... and label them with the respective list entries
               xticklabels=classes, yticklabels=classes,
               title=title,
               ylabel='True label',
               xlabel='Predicted label')

        # Rotate the tick labels and set their alignment.
        plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
                 rotation_mode="anchor")

        # Loop over data dimensions and create text annotations.
        fmt = '.2f' if normalize else 'd'
        thresh = cm.max() / 2.
        for i in range(cm.shape[0]):
            for j in range(cm.shape[1]):
                ax.text(j, i, format(cm[i, j], fmt),
                        ha="center", va="center",
                        color="white" if cm[i, j] > thresh else "black")
        fig.tight_layout()
        return ax

    plot_confusion_matrix(is_technical_gt, is_technical_infer_zeros)
    plt.savefig('/home/pierre/confusion.png')
    plt.show()


if __name__ == '__main__':
    test_zeros_classif()
