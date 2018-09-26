use_cuda = True
from scvi.dataset.dataset10X import Dataset10X
from scvi.dataset.pbmc import PbmcDataset
from scvi.harmonization.benchmark import assign_label
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent, run_model
import sys
from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer

model_type = str(sys.argv[1])
# option = str(sys.argv[2])
# plotname = 'Easy1'+option
plotname = 'Easy1'
print(model_type)

pbmc = PbmcDataset()
pbmc68k = Dataset10X('fresh_68k_pbmc_donor_a')

pbmc68k.cell_types = ['unlabelled']
pbmc68k.labels = np.repeat(0,len(pbmc68k)).reshape(len(pbmc68k),1)
gene_dataset = GeneExpressionDataset.concat_datasets(pbmc,pbmc68k)
gene_dataset.subsample_genes(5000)
print(gene_dataset.cell_types)

vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')
trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
trainer.train(n_epochs=250)

scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2)
scanvi.load_state_dict(vae.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)

trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices != 2))
trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
    indices=(gene_dataset.batch_indices == 2)
)
trainer_scanvi.train(n_epochs=50)

keys = gene_dataset.cell_types
latent, batch_indices, labels = trainer_scanvi.labelled_set.get_latent()
pred = trainer_scanvi.labelled_set.compute_predictions()
np.mean(pred[0] == pred[1])

scanvi_posterior = trainer_scanvi.create_posterior(trainer_scanvi.model,gene_dataset)
pred = scanvi_posterior.compute_predictions()

full = trainer.create_posterior(vae, gene_dataset, indices=np.arange(len(gene_dataset)))
scale1 = full.sequential().get_harmonized_scale(0)
scale2 = full.sequential().get_harmonized_scale(1)
np.save('../DE/de.scale1.npy',scale1)
np.save('../DE/de.scale2.npy',scale2)
np.save('../DE/labels.npy',full.gene_dataset.labels.ravel())
# full.gene_dataset.cell_types
# ['NK cells', 'B cells', 'CD4 T cells', 'unlabelled', 'Other',
#        'CD8 T cells', 'Megakaryocytes', 'Dendritic Cells',
#        'CD14+ Monocytes', 'FCGR3A+ Monocytes']
from scvi.inference.posterior import get_bayes_factors
bayes1 = get_bayes_factors(scale1,full.gene_dataset.labels.ravel(),0,4)
bayes2 = get_bayes_factors(scale2,full.gene_dataset.labels.ravel(),0,4)
gene_dataset.gene_names
