use_cuda = True
from scvi.harmonization.utils_chenling import get_matrix_from_dir
import numpy as np
from scvi.dataset.dataset import GeneExpressionDataset
from scvi.harmonization.utils_chenling import eval_latent
from scvi.models.vae import VAE
from scvi.inference import UnsupervisedTrainer
from copy import deepcopy


class MSDataset(GeneExpressionDataset):
    def __init__(self, patient='MS19270',celltype='PBMCs'):
        self.patient = patient
        self.celltype = celltype
        count, labels, cell_type, gene_names = self.preprocess()
        super(MSDataset, self).__init__(
            *GeneExpressionDataset.get_attributes_from_matrix(
                count, labels=labels),
            gene_names=np.char.upper(gene_names), cell_types=cell_type)
    def preprocess(self):
        print("Preprocessing data")
        count, geneid, cellid = get_matrix_from_dir('Gerd/{0}/{1}/GRCh38'.format(self.patient,self.celltype))
        geneid = geneid[:, 1]
        count = count.T.tocsr()
        print("Finish preprocessing data")
        labels = np.repeat(0,len(cellid))
        cell_type = ['unlabelled']
        return count,labels,cell_type,geneid

patient_id = ['MS19270', 'MS49131'  ,  'PTC32190',  'PTC41540',  'PTC85037']
dataset = []
for id in patient_id:
    data = MSDataset(id,'PBMCs')
    data.subsample_genes(data.nb_genes)
    dataset.append(data)
    data = MSDataset(id, 'CSF')
    data.subsample_genes(data.nb_genes)
    dataset.append(data)

gene_dataset = GeneExpressionDataset.concat_datasets(dataset[0],dataset[1],
                                                     dataset[2],dataset[3],
                                                     dataset[4],dataset[5],
                                                     dataset[6],dataset[7],
                                                     dataset[8],dataset[9])
#
# gene_dataset.X = sparse.csr_matrix(gene_dataset.X)
#
gene_dataset.subsample_genes(5000)
# gene_dataset.X = gene_dataset.X.todense()
vae = VAE(gene_dataset.nb_genes, n_batch=gene_dataset.n_batches, n_labels=gene_dataset.n_labels,
          n_hidden=128, n_latent=10, n_layers=2, dispersion='gene')
trainer = UnsupervisedTrainer(vae, gene_dataset, train_size=1.0)
trainer.train(n_epochs=100)

# if os.path.isfile('../' + filename + '/' + 'vae' + '.rep' + str(rep) + '.pkl'):
#     trainer.model = torch.load('../' + filename + '/' + 'vae' + '.rep' + str(rep) + '.pkl')
# else:
#     if gene_dataset.X.shape[0] > 40000:
#         trainer.train(n_epochs=100)
#     elif gene_dataset.X.shape[0] < 20000:
#         trainer.train(n_epochs=500)
#     else:
#         trainer.train(n_epochs=250)
#
import torch
torch.save(trainer.model, '../' + 'Gerd' + '/' + 'vae.pkl')

M_sampling=1
px_scales = []
all_labels = []

imputed_list=[]
import torch

trainer.train_set.get_harmonized_scale()

for tensors in trainer.train_set:
    sample_batch, local_l_mean, local_l_var, batch_index, label = tensors
    if trainer.model.log_variational:
        sample_batch = torch.log(1 + sample_batch)
    qz_m, qz_v, z = trainer.model.z_encoder(sample_batch)
    batch_index = torch.cuda.IntTensor(len(label),1).fill_(0)
    library = torch.cuda.FloatTensor(len(label),1).fill_(4)
    px_scale, px_r, px_rate, px_dropout = trainer.model.decoder('gene', z, library, batch_index)

from scvi.inference.posterior import get_bayes_factors
bayes_factors_list = get_bayes_factors(px_scale, all_labels, cell_idx, other_cell_idx=other_cell_idx,
                                               M_permutation=M_permutation, permutation=permutation)
from torch.distributions import Normal, kl_divergence as kl

