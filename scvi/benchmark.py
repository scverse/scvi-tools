from scvi.dataset import CortexDataset
from scvi.inference import VariationalInference
from scvi.models import VAE


def cortex_benchmark(n_epochs=250, use_cuda=True, unit_test=False):
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, use_cuda=use_cuda)
    infer_cortex_vae.train(n_epochs=n_epochs)

    infer_cortex_vae.ll('test')  # assert ~ 1200
    infer_cortex_vae.differential_expression('test')
    infer_cortex_vae.imputation('test', rate=0.1)  # assert ~ 2.3
    n_samples = 1000 if not unit_test else 10
    infer_cortex_vae.show_t_sne('test', n_samples=n_samples)
    return infer_cortex_vae


def benchmark(dataset, n_epochs=250, use_cuda=True):
    vae = VAE(dataset.nb_genes, n_batch=dataset.n_batches)
    infer = VariationalInference(vae, dataset, use_cuda=use_cuda)
    infer.train(n_epochs=n_epochs)
    infer.ll('test')
    infer.imputation('test', rate=0.1)  # assert ~ 2.1
    return infer


def harmonization_benchmarks(n_epochs=1, use_cuda=True):
    # retina_benchmark(n_epochs=n_epochs)
    pass


def annotation_benchmarks(n_epochs=1, use_cuda=True):
    # some cortex annotation benchmark
    pass


def all_benchmarks(n_epochs=250, use_cuda=True, unit_test=False):
    cortex_benchmark(n_epochs=n_epochs, use_cuda=use_cuda, unit_test=unit_test)

    harmonization_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
    annotation_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
