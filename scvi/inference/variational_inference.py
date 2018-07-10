import torch
from torch.nn import functional as F

from scvi.dataset.utils import TrainTestDataLoaders, AlternateSemiSupervisedDataLoaders, JointSemiSupervisedDataLoaders
from . import Inference, ClassifierInference


class VariationalInference(Inference):
    metrics = ['ll', 'accuracy', 'imputation', 'differential_expression', 'batch_entropy_mixing']
    tasks = ['imputation_stats', 'differential_expression_stats', 'show_t_sne']
    default_metrics_to_monitor = ['ll']

    def __init__(self, model, gene_dataset, train_size=0.1, **kwargs):
        super(VariationalInference, self).__init__(model, gene_dataset, **kwargs)
        self.kl = None
        self.data_loaders = TrainTestDataLoaders(self.gene_dataset, train_size=train_size, pin_memory=self.use_cuda)

    def loss(self, tensors):
        sample_batch, local_l_mean, local_l_var, batch_index, _ = tensors
        reconst_loss, kl_divergence = self.model(sample_batch, local_l_mean, local_l_var, batch_index)
        loss = torch.mean(reconst_loss + self.kl_weight * kl_divergence)
        return loss

    def on_epoch_begin(self):
        self.kl_weight = self.kl if self.kl is not None else min(1, self.epoch / self.n_epochs)


class SemiSupervisedVariationalInference(VariationalInference):
    metrics = VariationalInference.metrics + ['accuracy']
    default_metrics_to_monitor = VariationalInference.default_metrics_to_monitor + ['accuracy']
    baselines = ['svc_rf']


class AlternateSemiSupervisedVariationalInference(SemiSupervisedVariationalInference):
    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, n_epochs_classifier=1,
                 lr_classification=0.1, **kwargs):
        super(AlternateSemiSupervisedVariationalInference, self).__init__(model, gene_dataset, **kwargs)

        self.n_epochs_classifier = n_epochs_classifier
        self.lr_classification = lr_classification
        self.data_loaders = AlternateSemiSupervisedDataLoaders(gene_dataset, n_labelled_samples_per_class)

        self.classifier_inference = ClassifierInference(
            model.classifier, gene_dataset, metrics_to_monitor=[], benchmark=True,
            data_loaders=self.data_loaders.classifier_data_loaders(), sampling_model=self.model
        )

    def on_epoch_end(self):
        self.classifier_inference.fit(self.n_epochs_classifier, lr=self.lr_classification)
        return super(AlternateSemiSupervisedVariationalInference, self).on_epoch_end()


class JointSemiSupervisedVariationalInference(SemiSupervisedVariationalInference):
    def __init__(self, model, gene_dataset, n_labelled_samples_per_class=50, classification_ratio=100, **kwargs):
        super(JointSemiSupervisedVariationalInference, self).__init__(model, gene_dataset, **kwargs)
        self.data_loaders = JointSemiSupervisedDataLoaders(gene_dataset, n_labelled_samples_per_class)
        self.classification_ratio = classification_ratio

    def loss(self, tensors_unlabelled, tensors_labelled):
        loss = super(JointSemiSupervisedVariationalInference, self).loss(tensors_unlabelled)
        sample_batch, _, _, _, y = tensors_labelled
        classification_loss = F.cross_entropy(self.model.classify(sample_batch), y.view(-1))
        loss += classification_loss * self.classification_ratio
        return loss
