import types

import torch
import torch.nn.functional as F

from scvi.models.classifier import Classifier


def adversarial_loss(self, tensors, *next_tensors):
    if self.epoch > self.warm_up:
        sample_batch, _, _, batch_index, label = tensors
        qm_z, _, _ = self.model.z_encoder(torch.log(1 + sample_batch), label)  # label only used in VAEC
        cls_loss = (self.scale * F.cross_entropy(self.adversarial_cls(qm_z), batch_index.view(-1)))
        self.optimizer_cls.zero_grad()
        cls_loss.backward(retain_graph=True)
        self.optimizer_cls.step()
    else:
        cls_loss = 0
    return type(self).loss(self, tensors, *next_tensors) - cls_loss


def adversarial_loss_fish(self, tensors_seq, tensors_fish):
    if self.epoch > self.warm_up:
        sample_batch, local_l_mean, local_l_var, batch_index, labels = tensors_seq
        z = self.model.sample_from_posterior_z(sample_batch, mode="scRNA")
        cls_loss = (self.scale * F.cross_entropy(self.adversarial_cls(z), torch.zeros_like(batch_index).view(-1)))
        sample_batch_fish, local_l_mean, local_l_var, batch_index_fish, _, _, _ = tensors_fish
        z = self.model.sample_from_posterior_z(sample_batch, mode="smFISH")
        cls_loss += (self.scale * F.cross_entropy(self.adversarial_cls(z), torch.ones_like(batch_index).view(-1)))
        self.optimizer_cls.zero_grad()
        cls_loss.backward(retain_graph=True)
        self.optimizer_cls.step()
    else:
        cls_loss = 0
    return type(self).loss(self, tensors_seq, tensors_fish) - cls_loss


def adversarial_train(self, n_epochs=20, lr=1e-3, weight_decay=1e-4):
    self.adversarial_cls = Classifier(self.model.n_latent, n_labels=self.model.n_batch, n_layers=3)
    if self.use_cuda:
        self.adversarial_cls.cuda()
    self.optimizer_cls = torch.optim.Adam(filter(lambda p: p.requires_grad, self.adversarial_cls.parameters()), lr=lr,
                                          weight_decay=weight_decay)
    type(self).train(self, n_epochs=n_epochs, lr=lr)


def adversarial_wrapper(trainer, warm_up=100, scale=50, mode="regular"):
    trainer.warm_up = warm_up
    trainer.scale = scale
    if mode == "regular":
        trainer.loss = types.MethodType(adversarial_loss, trainer)
    else:
        trainer.loss = types.MethodType(adversarial_loss_fish, trainer)
    trainer.train = types.MethodType(adversarial_train, trainer)
    return trainer
