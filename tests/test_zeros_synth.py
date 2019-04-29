import numpy as np
import torch

from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
from scvi.dataset.synthetic import ZISyntheticDatasetCorr, SyntheticDatasetCorr, ZISyntheticDatasetCorrDistinct
from sklearn.metrics import confusion_matrix

from zifa_full import VAE as ZifaVae

import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def test_enough_zeros():
    """
    Can be seen as a 'pre' test for test_zeros_classif
    In test_zeros_classif, we classify the zeros obtained by scVI
    Hence, the objective of test_enough_zeros is to check that the synthetic
    dataset has enough zeros and that the proportion of technical zeros over biological zeros
    is somehow balanced
    :return:
    """

    # nb_data = ZISyntheticDatasetCorr(lam_0=180, n_cells_cluster=2000,
    #                                  weight_high=6, weight_low=3.5, n_genes_high=35,
    #                                  n_genes_total=50)
    nb_data = ZISyntheticDatasetCorrDistinct(lam_0=180, n_cells_cluster=2000,
                                                weight_high=6, weight_low=3.5, n_genes_high=35,
                                                n_genes_total=50,
                                                lam_dropout_low=0.0, lam_dropout_high=0.0)
    is_technical_mask = nb_data.is_technical.squeeze()
    nb_data_zeros = nb_data.X == 0
    tech_zeros = nb_data_zeros[is_technical_mask].sum()
    bio_zeros = nb_data_zeros[~is_technical_mask].sum()
    print("Prop of technical zeros :", tech_zeros)
    print("Prop of biological zeros :", bio_zeros)
    assert .2 <= tech_zeros / float(bio_zeros) <= 5.
    assert tech_zeros >= 1000


def test_zeros_classif():
    """
    Test that controls that scVI inferred distributions make sense on a non-trivial synthetic
    dataset.

    We define technical zeros of the synthetic dataset as the zeros that result from
    highly expressed genes (relatively to the considered cell) and the biological zeros as the
    rest of the zeros
    :return: None
    """
    torch.manual_seed(seed=42)
    # synth_data = ZISyntheticDatasetCorr(n_cells_cluster=200, n_clusters=3,
    #                                     n_genes_high=25, n_genes_total=50,
    #                                     weight_high=4e-2, weight_low=1e-2,
    #                                     lam_0=50., n_batches=1)
    # Step 1: Generating dataset and figuring out ground truth values
    synth_data = ZISyntheticDatasetCorrDistinct(lam_0=180, n_cells_cluster=2000,
                                                weight_high=6, weight_low=3.5, n_genes_high=15,
                                                n_genes_total=45)



    zeros_mask = (synth_data.X == 0).astype(bool)
    poisson_params_gt = synth_data.exprs_param.squeeze()

    is_technical_all = synth_data.is_technical[0, :]
    is_technical_gt = is_technical_all[zeros_mask]  # 1d array

    # Step 2: Training scVI model
    mdl = VAE(n_input=synth_data.nb_genes, n_batch=synth_data.n_batches,
              reconstruction_loss='zinb', n_latent=1)
    # mdl = ZifaVae(n_genes=synth_data.nb_genes, n_batch=synth_data.n_batches,
    #               reconstruction_loss='zinb')

    trainer = UnsupervisedTrainer(model=mdl, gene_dataset=synth_data, use_cuda=True, train_size=1.0,
                                  # frequency=1,
                                  # early_stopping_kwargs={
                                  #    'early_stopping_metric': 'll',
                                  #    'save_best_state_metric': 'll',
                                  #    'patience': 15,
                                  #    'threshold': 3}
                                  )
    trainer.train(n_epochs=250, lr=1e-3)
    full = trainer.create_posterior(trainer.model, synth_data,
                                    indices=np.arange(len(synth_data)))

    # Step 3: Inference
    technical_scores_infer = []
    poisson_params = []
    p_dropout_infered = []
    with torch.no_grad():
        for tensors in full.sequential():
            sample_batch, _, _, batch_index, labels = tensors
            px_scale, px_dispersion, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v, library = mdl.inference(
                sample_batch, batch_index)
            p_zero = 1.0 / (1.0 + torch.exp(-px_dropout))
            p_dropout_infered.append(p_zero.cpu().numpy())

            technical_counts_batch = torch.zeros((sample_batch.size(0), sample_batch.size(1), 100))
            zero_counts_batch = torch.zeros_like(technical_counts_batch)
            l_train_batch = torch.zeros_like(technical_counts_batch)

            for n_mc_sim in range(100):
                p = px_rate / (px_rate + px_dispersion)
                r = px_dispersion
                l_train = torch.distributions.Gamma(concentration=r, rate=(1 - p) / p).sample()
                l_train = torch.clamp(l_train, max=1e18)
                X = torch.distributions.Poisson(l_train).sample()
                l_train_batch[:, :, n_mc_sim] = l_train
                p_zero = 1.0 / (1.0 + torch.exp(-px_dropout))
                random_prob = torch.rand_like(p_zero)
                X[random_prob <= p_zero] = 0

                zero_counts_batch[:, :, n_mc_sim] = X == 0
                technical_counts_batch[:, :, n_mc_sim] = (random_prob <= p_zero)

            l_train_batch = torch.mean(l_train_batch, dim=(-1))
            poisson_params.append(l_train_batch)

            zero_counts_batch = torch.mean(zero_counts_batch, dim=(-1))
            technical_counts_batch = torch.mean(technical_counts_batch, dim=(-1))
            technical_score_batch = (technical_counts_batch / zero_counts_batch)
            technical_scores_infer.append(technical_score_batch.cpu().numpy())

    technical_scores_infer = np.concatenate(technical_scores_infer)
    technical_scores_infer[np.isnan(technical_scores_infer)] = 0.0
    technical_scores_zeros = technical_scores_infer[zeros_mask]
    is_technical_infer_zeros = technical_scores_infer[zeros_mask] >= 0.5
    assert technical_scores_infer.shape == zeros_mask.shape
    assert is_technical_infer_zeros.shape == is_technical_gt.shape

    # Final Step: Checking predictions
    # Technical pos
    fig, axes = plt.subplots(ncols=2, figsize=(10, 6))
    vmin = min(technical_scores_infer.min(), is_technical_all.min())
    vmax = max(technical_scores_infer.max(), is_technical_all.max())

    sns.heatmap(technical_scores_infer, vmin=vmin, vmax=vmax, ax=axes[1])
    axes[1].set_title('Tech Zeros Predicted')
    sns.heatmap(is_technical_all, vmin=vmin, vmax=vmax, ax=axes[0])
    axes[0].set_title('Tech Zeros GT')
    plt.savefig('tech_zeros.png')
    plt.close()

    sorted_indices = np.argsort(is_technical_gt)
    plt.plot(technical_scores_zeros[sorted_indices], label='Inferred technical counts')
    plt.plot(is_technical_gt[sorted_indices], label='Is Technical')
    plt.legend()
    plt.savefig('zeros_pos_dropout_corr.png')

    # Dropout checks
    p_dropout_infered_all = np.concatenate(p_dropout_infered)
    p_dropout_gt = synth_data.p_dropout.squeeze()
    # vmin = min(p_dropout_gt.min(), p_dropout_infered_all.min())
    # vmax = max(p_dropout_gt.max(), p_dropout_infered_all.max())
    vmin = 0.0
    vmax = 2.0*p_dropout_gt.max()
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
    sns.heatmap(p_dropout_infered_all, vmin=vmin, vmax=vmax, ax=axes[0, 1])
    axes[0, 1].set_title('Dropout Rate Predicted')

    sns.heatmap(p_dropout_gt, vmin=vmin, vmax=vmax, ax=axes[0, 0])
    axes[0, 0].set_title('Dropout Rate GT')

    # Poisson Params checks
    poisson_params = np.concatenate(poisson_params)
    vmin = min(poisson_params_gt.min(), poisson_params.min())
    vmax = max(poisson_params_gt.max(), poisson_params.max())
    # vmin = 0.0
    # vmax = 10.0
    sns.heatmap(poisson_params, vmin=vmin, vmax=vmax, ax=axes[1, 1])
    axes[1, 1].set_title('Poisson Distribution Parameter Predicted')

    sns.heatmap(poisson_params_gt,  vmin=vmin, vmax=vmax, ax=axes[1, 0])
    axes[1, 0].set_title('Poisson Distribution Parameter GT')
    plt.savefig('params_comparison.png')
    plt.close()

    print('Average Poisson L1 error: ', np.abs(poisson_params - poisson_params_gt).mean())
    sns.heatmap(np.abs(poisson_params - poisson_params_gt) / poisson_params_gt)
    plt.savefig('params_diff.png')
    plt.close()

    # Tech/Bio Classif checks
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
        return ax, cm
    ax, _ = plot_confusion_matrix(is_technical_gt, is_technical_infer_zeros, normalize=False)
    plt.savefig('confusion.png')
    plt.close()

    ax, cm = plot_confusion_matrix(is_technical_gt, is_technical_infer_zeros, normalize=True)
    plt.savefig('confusion_normalized.png')
    plt.close()

    assert cm[0, 0] >= 0.6
    assert cm[1, 1] >= 0.6
    assert cm[1, 0] <= 0.4
    assert cm[0, 1] <= 0.4


if __name__ == '__main__':
    test_zeros_classif()
