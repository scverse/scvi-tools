import scvi
from scvi.core.models.vae import VAE
from scvi.core.trainers.inference import UnsupervisedTrainer
from scvi.core.trainers.annotation import ClassifierTrainer
from scvi.core.models.classifier import Classifier

scvi.set_seed(0)
use_cuda = True


def test_sampling_zl(save_path):
    cortex_dataset = scvi.dataset.cortex(save_path=save_path)
    scvi.dataset.setup_anndata(cortex_dataset, labels_key="cell_type")
    cortex_vae = VAE(
        cortex_dataset.uns["_scvi"]["summary_stats"]["n_genes"],
        cortex_dataset.uns["_scvi"]["summary_stats"]["n_batch"],
    )
    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae, cortex_dataset, train_size=0.5, use_cuda=use_cuda
    )
    trainer_cortex_vae.train(n_epochs=2)

    cortex_cls = Classifier(
        (cortex_vae.n_latent + 1),
        n_labels=cortex_dataset.uns["_scvi"]["summary_stats"]["n_labels"],
    )
    trainer_cortex_cls = ClassifierTrainer(
        cortex_cls, cortex_dataset, sampling_model=cortex_vae, sampling_zl=True
    )
    trainer_cortex_cls.train(n_epochs=2)
    trainer_cortex_cls.test_set.accuracy()


def test_annealing_procedures(save_path):
    cortex_dataset = scvi.dataset.cortex(save_path=save_path)
    scvi.dataset.setup_anndata(cortex_dataset, labels_key="cell_type")

    cortex_vae = VAE(
        cortex_dataset.uns["_scvi"]["summary_stats"]["n_genes"],
        cortex_dataset.uns["_scvi"]["summary_stats"]["n_batch"],
    )

    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae,
        cortex_dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        n_epochs_kl_warmup=1,
    )
    trainer_cortex_vae.train(n_epochs=2)
    assert trainer_cortex_vae.kl_weight >= 0.99, "Annealing should be over"

    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae,
        cortex_dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        n_epochs_kl_warmup=5,
    )
    trainer_cortex_vae.train(n_epochs=2)
    assert trainer_cortex_vae.kl_weight <= 0.99, "Annealing should be proceeding"

    # iter
    trainer_cortex_vae = UnsupervisedTrainer(
        cortex_vae,
        cortex_dataset,
        train_size=0.5,
        use_cuda=use_cuda,
        n_iter_kl_warmup=1,
        n_epochs_kl_warmup=None,
    )
    trainer_cortex_vae.train(n_epochs=2)
    assert trainer_cortex_vae.kl_weight >= 0.99, "Annealing should be over"
