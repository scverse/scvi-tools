use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir
from scvi.harmonization.benchmark import assign_label
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.dataset.pbmc import Dataset10X
from scvi.dataset.cite_seq import CiteSeqDataset

from scvi.harmonization.utils_chenling import eval_latent, run_model
from scvi.inference.posterior import *
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
use_cuda = True
from scvi.metrics.clustering import select_indices_evenly, clustering_scores, entropy_batch_mixing
from scipy import sparse
from scvi.models.vae import VAE
from scvi.models.scanvi import SCANVI
from scvi.inference import UnsupervisedTrainer, SemiSupervisedTrainer

cite = CiteSeqDataset('pbmc')
cite.subsample_genes(cite.nb_genes)
cite.gene_names = cite.gene_symbols
cite.cell_types = ['unlabelled']

donner = Dataset10X('fresh_68k_pbmc_donor_a')
donner.cell_types = np.asarray(['unlabelled'])
donner.subsample_genes(donner.nb_genes)
donner.gene_names = donner.gene_symbols

cell_types = np.array(["cd4_t_helper", "regulatory_t", "naive_t", "memory_t", "cytotoxic_t", "naive_cytotoxic",
                       "b_cells", "cd4_t_helper", "cd34", "cd56_nk", "cd14_monocytes"])
datasets = []
for cell_type in cell_types:
    dataset = Dataset10X(cell_type, save_path='data/')
    dataset.cell_types = np.array([cell_type])
    dataset.subsample_genes(dataset.nb_genes)
    dataset.gene_names = dataset.gene_symbols
    datasets += [dataset]

pure = GeneExpressionDataset.concat_datasets(*datasets, shared_batches=True)


plotname = 'PureDonorCite.vae'
gene_dataset = GeneExpressionDataset.concat_datasets(cite, donner, pure)
gene_dataset.X = sparse.csr_matrix(gene_dataset.X)
gene_dataset.subsample_genes(5000)
gene_dataset.X = gene_dataset.X.todense()


vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')
trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
trainer.train(n_epochs=250)
batch_entropy = trainer.train_set.entropy_batch_mixing()
keys = gene_dataset.cell_types
vae_full = trainer.create_posterior(trainer.model, gene_dataset, indices=np.arange(len(gene_dataset)))
vae_latent, batch_indices, labels = vae_full.sequential().get_latent()
batch_indices = batch_indices.ravel()
labels = labels.ravel()


scanvi = SCANVI(gene_dataset.nb_genes, gene_dataset.n_batches, gene_dataset.n_labels, n_layers=2)
scanvi.load_state_dict(vae.state_dict(), strict=False)
trainer_scanvi = SemiSupervisedTrainer(scanvi, gene_dataset, classification_ratio=1,
                                       n_epochs_classifier=1, lr_classification=5 * 1e-3)

trainer_scanvi.labelled_set = trainer_scanvi.create_posterior(indices=(gene_dataset.batch_indices == 2))
trainer_scanvi.unlabelled_set = trainer_scanvi.create_posterior(
    indices=(gene_dataset.batch_indices != 2)
)
trainer_scanvi.train(n_epochs=50)
scanvi_posterior = trainer_scanvi.create_posterior(trainer_scanvi.model,gene_dataset)
pred = scanvi_posterior.sequential().compute_predictions()
np.mean(pred[0][batch_indices==2] == pred[1][batch_indices==2])

keys = gene_dataset.cell_types
latent, batch_indices, labels = trainer_scanvi.unlabelled_set.get_latent()
batch_indices = batch_indices.ravel()


scale1 = full.sequential().get_harmonized_scale(0)
scale2 = full.sequential().get_harmonized_scale(1)
from scvi.inference.posterior import get_bayes_factors
bayes1 = get_bayes_factors(scale1,pred[1],0,5)
bayes2 = get_bayes_factors(scale2,pred[1],0,5)
np.savetxt('../genenames.txt',gene_dataset.gene_names,fmt='%s',delimiter=',')
np.savetxt('../keys.txt',keys,fmt='%s',delimiter=',')

np.save('../three2.latent.npy',latent)
np.save('../three2.batch_indices.npy',batch_indices)
np.save('../three2.labels.npy',labels)
np.savetxt('../three2.keys.txt',keys,fmt='%s',delimiter=',')

latent = np.load('../three/three2.latent.npy')
batch_indices = np.load('../three/three2.batch_indices.npy')
labels = np.load('../three/three2.labels.npy')
keys = np.genfromtxt('../three/three2.keys.txt',dtype=str)

# 27.52s/it
sample = select_indices_evenly(2000, batch_indices)
res = clustering_scores(np.asarray(latent)[sample, :], labels[sample], 'knn', len(np.unique(labels[sample])))
colors = sns.color_palette('tab20')
latent_s = latent[sample, :]
label_s = labels[sample]
batch_s = batch_indices[sample]
if latent_s.shape[1] != 2:
    latent_s = TSNE().fit_transform(latent_s)


plt.figure(figsize=(10, 10))
plt.scatter(latent_s[:, 0], latent_s[:, 1], c=batch_s, edgecolors='none')
plt.axis("off")
plt.tight_layout()
plt.savefig('../' + plotname + '.batchid.png')

def rank(x):
    uniq = np.unique(x)
    lookup = dict(zip(uniq,np.arange(len(uniq))))
    x = np.asarray([lookup[i] for i in x])
    return(x)

for batch in [0,1,2]:
    latent_sb = latent_s[batch_s==batch,:]
    label_sb = label_s[batch_s==batch]
    fig, ax = plt.subplots(figsize=(13, 10))
    key_sb = keys[np.unique(label_sb).astype('int')]
    label_sb = rank(label_sb)
    key_order = np.argsort(key_sb)
    for i, k in enumerate(key_order):
        ax.scatter(latent_sb[label_sb == k, 0], latent_sb[label_sb == k, 1], c=colors[i % 20], label=key_sb[k],
                   edgecolors='none')
        ax.legend(bbox_to_anchor=(1.1, 0.5), borderaxespad=0, fontsize='x-large')
    fig.tight_layout()
    fig.savefig('../' + plotname +'.'+ str(batch)+ '.labels.png')

sample = select_indices_evenly(2000, batch_indices)
batch_entropy = entropy_batch_mixing(latent[sample, :], batch_indices[sample])
print("Entropy batch mixing :", batch_entropy)

from sklearn.neighbors import NearestNeighbors
# 10 nearest neighbor in citeseq from each cell in pure population
nbrs = NearestNeighbors(n_neighbors=3 + 1).fit(latent[batch_indices==0,:])
indices = nbrs.kneighbors(latent[batch_indices==2,:], return_distance=False)[:, 1:]

start_freq = np.unique(labels[batch_indices==0],return_counts=True)
start_freq = dict(zip(keys[start_freq[0]], start_freq[1]/np.sum(start_freq[1])))


enrichment=[]
for j in np.unique(labels[batch_indices==2]):
    idx = indices[labels[batch_indices==2] == j]
    matched_lab = labels[batch_indices==0][idx.ravel()]
    freq = np.unique(matched_lab,return_counts=True)
    matched_names = keys[freq[0]]
    freq = freq[1]/np.sum(freq[1])
    res = dict(zip(matched_names,[freq[i]/start_freq[x] for i,x in enumerate(matched_names)]))
    enrichment.append([0 if x not in res.keys() else res[x] for i,x in enumerate(start_freq.keys())])

enrichment =np.array(enrichment)
np.argmax(enrichment[0])
np.apply_along_axis(np.argmax,1,enrichment)


adt_expression_clr = np.genfromtxt('../three/PBMC.ADT_cut_clr.txt').T
adt_expression_clr[adt_expression_clr<0]=0
from sklearn.preprocessing import scale
adt_expression_clr = scale(adt_expression_clr,with_mean=False)

markers = ["CD3","CD4","CD8","CD2","CD45RA","CD57","CD16","CD14","CD11c","CD19"]


protein=[]
for j in np.unique(labels[batch_indices==2]):
    idx = indices[labels[batch_indices==2] == j]
    clr = adt_expression_clr[idx.ravel()]
    protein.append(np.mean(clr,axis=0))

protein = np.asarray(protein)


def Heatmap(matrix, rownames,colnames,title):
    fig, ax = plt.subplots(figsize=(7,7))
    # We want to show all ticks...
    im = ax.imshow(matrix,aspect='auto')
    ax.set_xticks(np.arange(len(rownames)))
    ax.set_yticks(np.arange(len(colnames)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(rownames)
    ax.set_yticklabels(colnames)
    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")
    # Loop over data dimensions and create text annotations.
    for i in range(len(colnames)):
        for j in range(len(rownames)):
            text = ax.text(j, i, "{:.2f}".format(matrix[i, j]),
                           ha="center", va="center", color="w")
    ax.set_title(title)
    fig.tight_layout()
    plt.savefig('../heatmap.png')


Heatmap(enrichment,list(start_freq.keys()), keys[np.unique(labels[batch_indices==2])],"enrichment of cell type frequency")

Heatmap(protein,markers, keys[np.unique(labels[batch_indices==2])],'Average Scaled Protein Expression')
