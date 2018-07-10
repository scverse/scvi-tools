from scvi.dataset import CortexDataset
from scvi.models import VAE
from scvi.inference import VariationalInference


def cortex_benchmark():
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, train_size=0.1)
    infer_cortex_vae.fit(n_epochs=250)

    infer_cortex_vae.ll('test')
    infer_cortex_vae.de('test')
    infer_cortex_vae.imputation('test', rate=0.1)
    infer_cortex_vae.show_t_sne('test', n_samples=1000)
    return infer_cortex_vae
