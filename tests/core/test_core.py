import scvi
from scvi.core.models.vae import VAE
from scvi.core.models.vaec import VAEC
from scvi.core.trainers import (
    UnsupervisedTrainer,
    ClassifierTrainer,
    SemiSupervisedTrainer,
)
from scvi.core.models.classifier import Classifier
from scvi.core.trainers.inference import AdapterTrainer

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


def test_adapter_trainer():

    n_latent = 5
    adata = scvi.dataset.synthetic_iid()
    model = scvi.models.SCVI(adata, n_latent=n_latent)
    model.train(1, train_size=0.5)

    trainer = AdapterTrainer(model.model, adata, model.trainer.test_set)
    trainer.train(n_epochs=1, n_path=1)


def test_classifier_accuracy(save_path):
    cortex_dataset = scvi.dataset.cortex(save_path=save_path)
    scvi.dataset.setup_anndata(cortex_dataset, labels_key="labels")
    cls = Classifier(
        cortex_dataset.uns["_scvi"]["summary_stats"]["n_genes"],
        n_labels=cortex_dataset.uns["_scvi"]["summary_stats"]["n_labels"],
    )
    cls_trainer = ClassifierTrainer(
        cls,
        cortex_dataset,
        metrics_to_monitor=["accuracy"],
        frequency=1,
        early_stopping_kwargs={
            "early_stopping_metric": "accuracy",
            "save_best_state_metric": "accuracy",
        },
    )
    cls_trainer.train(n_epochs=2)
    cls_trainer.train_set.accuracy()


def test_vaec():
    synthetic_dataset = scvi.dataset.synthetic_iid()
    scvi.dataset.setup_anndata(
        synthetic_dataset, batch_key="batch", labels_key="labels"
    )
    stats = synthetic_dataset.uns["_scvi"]["summary_stats"]

    vaec = VAEC(stats["n_genes"], stats["n_batch"], stats["n_labels"])
    trainer_synthetic_vaec = SemiSupervisedTrainer(
        vaec,
        synthetic_dataset,
        use_cuda=use_cuda,
        frequency=1,
        early_stopping_kwargs={
            "early_stopping_metric": "reconstruction_error",
            "on": "labelled_set",
            "save_best_state_metric": "reconstruction_error",
        },
    )
    trainer_synthetic_vaec.train(n_epochs=2)


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
