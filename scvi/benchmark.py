from scvi.dataset import CortexDataset
from scvi.models import VAE
from scvi.inference import VariationalInference


def cortex_benchmark():
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes, cortex_dataset.n_labels, n_batch=cortex_dataset.n_batches)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, train_size=0.75)
    infer_cortex_vae.fit(n_epochs=400)
    infer_cortex_vae.ll('test')
    infer_cortex_vae.de('test')
    infer_cortex_vae.imputation('train', rate=0.1)
    infer_cortex_vae.show_t_sne('test', n_samples=500)
    return infer_cortex_vae
