import numpy as np
import torch
import os

from scvi.models import VAE
from scvi.inference import UnsupervisedTrainer
from scvi.dataset.synthetic import ZISyntheticDatasetCorr
from sklearn.manifold import TSNE
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
    nb_data = ZISyntheticDatasetCorr(n_clusters=8, n_genes_high=15, n_overlap=3,
                                     lam_0=320, n_cells_cluster=100,
                                     weight_high=1.714286, weight_low=1,
                                     dropout_coef_high=0.08, dropout_coef_low=0.05)

    print(nb_data.X.shape)
    print(nb_data.exprs_param.min(), nb_data.exprs_param.max())
    is_technical_mask = nb_data.is_technical.squeeze()
    nb_data_zeros = nb_data.X == 0
    tech_zeros = nb_data_zeros[is_technical_mask].sum()
    bio_zeros = nb_data_zeros[~is_technical_mask].sum()
    print("Prop of technical zeros :", tech_zeros)
    print("Prop of biological zeros :", bio_zeros)
    assert .2 <= tech_zeros / float(bio_zeros) <= 5.
    assert tech_zeros >= 1000


def test_model_fit(model_fit: bool):
    """
    Test that controls that scVI inferred distributions make sense on a non-trivial synthetic
    dataset.

    We define technical zeros of the synthetic dataset as the zeros that result from
    highly expressed genes (relatively to the considered cell) and the biological zeros as the
    rest of the zeros
    :return: None
    """
    print('model_fit set to : ', model_fit)
    folder = '/tmp/scVI_zeros_test'
    print('Saving graphs in : {}'.format(folder))
    if not os.path.exists(folder):
        os.makedirs(folder)

    n_epochs = 150 if model_fit else 1
    n_mc_sim_total = 100 if model_fit else 1
    n_cells_cluster = 1000 if model_fit else 100

    torch.manual_seed(seed=42)
    synth_data = ZISyntheticDatasetCorr(n_clusters=8, n_genes_high=15, n_overlap=8,
                                        lam_0=320, n_cells_cluster=n_cells_cluster,
                                        weight_high=1.714286, weight_low=1,
                                        dropout_coef_low=0.08, dropout_coef_high=0.05)

    is_high = synth_data.is_highly_exp.squeeze()
    poisson_params_gt = synth_data.exprs_param.squeeze()

    # Step 2: Training scVI model
    mdl = VAE(n_input=synth_data.nb_genes, n_batch=synth_data.n_batches,
              reconstruction_loss='zinb', n_latent=5)

    trainer = UnsupervisedTrainer(model=mdl, gene_dataset=synth_data, use_cuda=True, train_size=1.0)
    trainer.train(n_epochs=n_epochs, lr=1e-3)
    full = trainer.create_posterior(trainer.model, synth_data,
                                    indices=np.arange(len(synth_data)))

    # Step 3: Inference
    poisson_params = []
    p_dropout_infered = []
    latent_reps = []
    bio_zero_p = []
    tech_zero_p = []
    with torch.no_grad():
        for tensors in full.sequential():
            # TODO: Properly sample posterior
            sample_batch, _, _, batch_index, labels = tensors
            px_scale, px_dispersion, px_rate, px_dropout, qz_m, qz_v, z, ql_m, ql_v, library = mdl.inference(
                sample_batch, batch_index)
            p_zero = 1.0 / (1.0 + torch.exp(-px_dropout))
            p_dropout_infered.append(p_zero.cpu().numpy())

            l_train_batch = torch.zeros((sample_batch.size(0), sample_batch.size(1), n_mc_sim_total),
                                        device=sample_batch.device)

            for n_mc_sim in range(n_mc_sim_total):
                p = px_rate / (px_rate + px_dispersion)
                r = px_dispersion
                l_train = torch.distributions.Gamma(concentration=r, rate=(1 - p) / p).sample()
                l_train = torch.clamp(l_train, max=1e18)
                X = torch.distributions.Poisson(l_train).sample()
                l_train_batch[:, :, n_mc_sim] = l_train
                p_zero = 1.0 / (1.0 + torch.exp(-px_dropout))
                random_prob = torch.rand_like(p_zero)
                X[random_prob <= p_zero] = 0

            l_train_batch = torch.mean(l_train_batch, dim=(-1))

            bio_zero_prob_batch = torch.exp(-l_train_batch)
            tech_zero_prob_batch = p_zero

            bio_zero_p.append(bio_zero_prob_batch.cpu().numpy())
            tech_zero_p.append(tech_zero_prob_batch.cpu().numpy())
            latent_reps.append(z.cpu().numpy())
            poisson_params.append(l_train_batch.cpu().numpy())

    latent_reps = np.concatenate(latent_reps)
    bio_zero_p = np.concatenate(bio_zero_p)
    tech_zero_p = np.concatenate(tech_zero_p)
    bio_zero_tech_no = bio_zero_p * (1.0 - tech_zero_p)
    tech_zero_bio_no = (1.0 - bio_zero_p) * tech_zero_p

    # Final Step: Checking predictions
    # Dropout checks
    p_dropout_infered_all = np.concatenate(p_dropout_infered)
    p_dropout_gt = synth_data.p_dropout.squeeze()
    vmin = 0.0
    vmax = 2.0 * p_dropout_gt.max()
    fig, axes = plt.subplots(ncols=2, nrows=2, figsize=(10, 10))
    sns.heatmap(p_dropout_infered_all, vmin=vmin, vmax=vmax, ax=axes[0, 1])
    axes[0, 1].set_title('Dropout Rate Predicted')
    sns.heatmap(p_dropout_gt, vmin=vmin, vmax=vmax, ax=axes[0, 0])
    axes[0, 0].set_title('Dropout Rate GT')

    # Poisson Params checks
    poisson_params = np.concatenate(poisson_params)
    vmin = min(poisson_params_gt.min(), poisson_params.min())
    vmax = max(poisson_params_gt.max(), poisson_params.max())
    sns.heatmap(poisson_params, vmin=vmin, vmax=vmax, ax=axes[1, 1])
    axes[1, 1].set_title('Poisson Distribution Parameter Predicted')

    sns.heatmap(poisson_params_gt, vmin=vmin, vmax=vmax, ax=axes[1, 0])
    axes[1, 0].set_title('Poisson Distribution Parameter GT')
    plt.savefig(os.path.join(folder, 'params_comparison.png'))
    plt.close()

    # TODO: Decrease test tolerances
    l1_poisson = np.abs(poisson_params - poisson_params_gt).mean()
    if model_fit:
        print('Average Poisson L1 error: ', l1_poisson)
        assert l1_poisson <= 0.75, \
            'High Error on Poisson parameter inference'
        l1_dropout = np.abs(p_dropout_infered_all - synth_data.p_dropout).mean()
        print('Average Dropout L1 error: ', l1_dropout)
        assert l1_dropout <= 5e-2, \
            'High Error on Dropout parameter inference'

    # tSNE plot
    print("Computing tSNE rep ...")
    x_rep = TSNE(n_components=2).fit_transform(latent_reps)
    print("Done!")
    pos = np.random.permutation(len(x_rep))[:1000]
    labels = ['c_{}'.format(idx) for idx in synth_data.labels[pos].squeeze()]
    sns.scatterplot(x=x_rep[pos, 0], y=x_rep[pos, 1], hue=labels,
                    palette='Set2')
    plt.title('Synthetic Dataset latent space')
    plt.savefig(os.path.join(folder, 't_sne.png'))
    plt.close()

    # Tech/Bio Classif checks
    # --For high expressed genes
    # ---Poisson nul and ZI non null
    print(bio_zero_tech_no[is_high].mean(), synth_data.probas_zero_bio_tech_high[1, 0])
    # ---Poisson non nul and .
    print(tech_zero_bio_no[is_high].mean(), synth_data.probas_zero_bio_tech_high[0, 1])

    # --Low expressed expressend
    # ---Poisson nul and ZI non null
    print(bio_zero_tech_no[~is_high].mean(), synth_data.probas_zero_bio_tech_low[1, 0])
    # ---Poisson non nul and .
    print(tech_zero_bio_no[~is_high].mean(), synth_data.probas_zero_bio_tech_low[0, 1])

    diff1 = np.abs(bio_zero_tech_no[is_high].mean() - synth_data.probas_zero_bio_tech_high[1, 0])
    diff2 = np.abs(tech_zero_bio_no[is_high].mean() - synth_data.probas_zero_bio_tech_high[0, 1])
    diff3 = np.abs(bio_zero_tech_no[~is_high].mean() - synth_data.probas_zero_bio_tech_low[1, 0])
    diff4 = np.abs(tech_zero_bio_no[~is_high].mean() - synth_data.probas_zero_bio_tech_low[0, 1])

    if model_fit:
        assert diff1 <= 2e-2
        assert diff2 <= 2e-2
        assert diff3 <= 2e-2
        assert diff4 <= 2e-2
