from scvi.dataset import CortexDataset, BrainLargeDataset, RetinaDataset
from scvi.inference import VariationalInference
from scvi.models import VAE


def cortex_benchmark(n_epochs=250, use_cuda=True):
    cortex_dataset = CortexDataset()
    vae = VAE(cortex_dataset.nb_genes)
    infer_cortex_vae = VariationalInference(vae, cortex_dataset, train_size=0.1, use_cuda=use_cuda)
    infer_cortex_vae.fit(n_epochs=n_epochs)

    infer_cortex_vae.ll('test')  # assert ~ 1200
    infer_cortex_vae.differential_expression('test')
    infer_cortex_vae.imputation('test', rate=0.1)  # assert ~ 2.3
    # binomial perturbation scheme Figure 11
    infer_cortex_vae.show_t_sne('test', n_samples=1000)
    return infer_cortex_vae


def brain_large_benchmark(n_epochs=250, use_cuda=True):
    brain_large_dataset = BrainLargeDataset()
    vae = VAE(brain_large_dataset.nb_genes)
    infer = VariationalInference(vae, brain_large_dataset, train_size=0.1, use_cuda=use_cuda)
    infer.fit(n_epochs=n_epochs)
    infer.ll('test')
    infer.imputation('test', rate=0.1)  # assert ~ 2.1
    return infer


def retina_benchmark(n_epochs=250, use_cuda=True):
    brain_large_dataset = RetinaDataset()
    vae = VAE(brain_large_dataset.nb_genes)
    infer = VariationalInference(vae, brain_large_dataset, train_size=0.1, use_cuda=use_cuda)
    infer.fit(n_epochs=n_epochs)
    infer.batch_entropy_mixing('test')  # Figure 8
    infer.imputation('test', rate=0.1)  # binomial perturbation scheme Fig 11
    return infer


def hemato_benchmark(n_epochs=250, use_cuda=True):
    brain_large_dataset = RetinaDataset()
    vae = VAE(brain_large_dataset.nb_genes)
    infer = VariationalInference(vae, brain_large_dataset, train_size=0.1, use_cuda=use_cuda)
    infer.fit(n_epochs=n_epochs)
    infer.batch_entropy_mixing('test')
    infer.imputation('test', rate=0.1)
    # Fig. 9(d) - uniform perturbation scheme
    # Fig 11 - binomial perturbation scheme
    return infer


def harmonization_benchmarks(n_epochs=1, use_cuda=True):
    # retina_benchmark(n_epochs=n_epochs)
    pass


def annotation_benchmarks(n_epochs=1, use_cuda=True):
    # some cortex annotation benchmark
    pass


def all_benchmarks(n_epochs=250, use_cuda=True):
    cortex_benchmark(n_epochs=n_epochs, use_cuda=use_cuda)
    hemato_benchmark(n_epochs=n_epochs, use_cuda=use_cuda)
    brain_large_benchmark(n_epochs=n_epochs, use_cuda=use_cuda)
    # retina_benchmark(n_epochs=n_epochs) the user should have the retina dataset

    harmonization_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
    annotation_benchmarks(n_epochs=n_epochs, use_cuda=use_cuda)
