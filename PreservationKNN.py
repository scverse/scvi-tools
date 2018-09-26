from scvi.harmonization.utils_chenling import run_model, eval_latent
import matplotlib.pyplot as plt
from scvi.dataset.dataset import GeneExpressionDataset
import numpy as np
from scvi.dataset.muris_tabula import TabulaMuris
dataset1 = TabulaMuris('facs')
dataset2 = TabulaMuris('droplet')

gene_dataset = GeneExpressionDataset.concat_datasets(dataset1, dataset2)



# dataset1.X = dataset1.X.tocsr()
# dataset2.X = dataset2.X.tocsr()
# seurat = SEURAT()
# seurat.create_seurat(dataset1, 1)
# seurat.create_seurat(dataset2, 2)
# latent1, latent2 = seurat.get_pcs()
# latent, batch_indices, labels, keys = run_model('Seurat', dataset1, dataset2, gene_dataset, filename='Tech2', ngenes=5000)


latent1_seurat,latent2_seurat, latent_seurat, _ = run_model('SeuratPC', gene_dataset, dataset1, dataset2, ngenes=5000)
KNeighbors  = np.concatenate([np.arange(10,100,10),np.arange(100,1000,100),np.arange(1000,2000,1000)])
batch_indices = gene_dataset.batch_indices.ravel()

# np.save('../cont/PBMC.latent1.seurat.npy', arr= latent1_seurat)
# np.save('../cont/PBMC.latent2.seurat.npy', arr= latent2_seurat)
# np.save('../cont/PBMC.latent.seurat.npy', arr= latent_seurat)
# latent1_seurat = np.load('../cont/TM.latent1.seurat.npy' )
# latent2_seurat = np.load('../cont/TM.latent2.seurat.npy')
# latent_seurat = np.load('../cont/TM.latent.seurat.npy')


latent1, _, label1, _ = run_model('vae', dataset1, 0, 0, ngenes=5000)
latent2, _, label2, _ = run_model('vae', dataset2, 0, 0, ngenes=5000)
latent, batch_indices, labels, keys = run_model('vae', gene_dataset, dataset1, dataset2,ngenes=5000)
eval_latent(batch_indices, labels, latent, keys)
KNeighbors  = np.concatenate([np.arange(10,100,10),np.arange(100,1000,100),np.arange(1000,2000,1000)])
np.save('../cont/PBMC.latent1.vae.npy', arr= latent1)
np.save('../cont/PBMC.latent2.vae.npy', arr= latent2)
np.save('../cont/PBMC.latent.vae.npy', arr= latent)

res_seurat = [KNNJaccardIndex(latent1_seurat,latent2_seurat,latent_seurat,batch_indices,i)[0] for i in KNeighbors]
res_vae = [KNNJaccardIndex(latent1,latent2,latent,batch_indices,i)[0] for i in KNeighbors]
res1 = [KNNJaccardIndex(latent1_seurat,latent2_seurat,latent,batch_indices,i)[0] for i in KNeighbors]
res2 = [KNNJaccardIndex(latent1,latent2,latent_seurat,batch_indices,i)[0] for i in KNeighbors]

plt.figure(figsize=(10, 10))
plt.plot(KNeighbors, res_seurat,'r',label='Seurat_Seurat')
plt.plot(KNeighbors, res_vae,'b',label='VAE_VAE')
plt.plot(KNeighbors, res1,'g',label='VAE_Seurat')
plt.plot(KNeighbors, res2,'y',label='Seurat_VAE')
legend = plt.legend(loc='lower right', shadow=False)
plt.savefig('../Easy1.KNN.png')
