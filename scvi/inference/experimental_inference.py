import types
from math import sqrt, pi

import torch
import torch.nn.functional as F
from torch.distributions import Normal, Uniform

from scvi.models.classifier import Classifier


def mmd_fourier(x1, x2, bandwidth=2., dim_r=500):
    d = x1.size(1)
    rw_n = sqrt(2. / bandwidth) * Normal(0., 1. / sqrt(d)).sample((dim_r, d)).type(x1.type())
    rb_u = 2 * pi * Uniform(0., 1.).sample((dim_r,)).type(x1.type())
    rf0 = sqrt(2. / dim_r) * torch.cos(F.linear(x1, rw_n, rb_u))
    rf1 = sqrt(2. / dim_r) * torch.cos(F.linear(x2, rw_n, rb_u))
    result = (torch.pow(rf0.mean(dim=0) - rf1.mean(dim=0), 2)).sum()
    return torch.sqrt(result)


def mmd_objective(z, batch_index, n_batch):
    mmd_method = mmd_fourier

    z_dim = z.size(1)
    batch_index = batch_index.view(-1)

    # STEP 1: construct lists of samples in their proper batches
    z_part = [z[batch_index == b_i] for b_i in range(n_batch)]

    # STEP 2: add noise to all of them and get the mmd
    mmd = 0
    for j, z_j in enumerate(z_part):
        z0_ = z_j
        aux_z0 = Normal(0., 1.).sample((1, z_dim)).type(z0_.type())
        z0 = torch.cat((z0_, aux_z0), dim=0)
        if len(z_part) == 2:
            z1_ = z_part[j + 1]
            aux_z1 = Normal(0., 1.).sample((1, z_dim)).type(z1_.type())
            z1 = torch.cat((z1_, aux_z1), dim=0)
            return mmd_method(z0, z1)
        z1 = z
        mmd += mmd_method(z0, z1)
    return mmd


def mmd_loss(self, tensors, *next_tensors):
    if self.epoch > self.warm_up:  # Leave a warm-up
        sample_batch, _, _, batch_index, label = tensors
        qm_z, _, _ = self.model.z_encoder(torch.log(1 + sample_batch), label)  # label only used in VAEC
        loss = mmd_objective(qm_z, batch_index, self.gene_dataset.n_batches)
    return type(self).loss(self, tensors, *next_tensors) + loss


def mmd_wrapper(infer, warm_up=100, scale=50):
    infer.warm_up = warm_up
    infer.scale = scale
    infer.loss = types.MethodType(adversarial_loss, infer)
    infer.train = types.MethodType(adversarial_train, infer)
    return infer


def adversarial_loss(self, tensors, *next_tensors):
    if self.epoch > self.warm_up:
        sample_batch, _, _, batch_index, label = tensors
        qm_z, _, _ = self.model.z_encoder(torch.log(1 + sample_batch), label)  # label only used in VAEC
        cls_loss = (self.scale * F.cross_entropy(self.GAN1(qm_z), batch_index.view(-1)))
        self.optimizer_GAN.zero_grad()
        cls_loss.backward(retain_graph=True)
        self.optimizer_GAN.step()
    else:
        cls_loss = 0
    return type(self).loss(self, tensors, *next_tensors) - cls_loss


def adversarial_train(self, n_epochs=20, lr=1e-3, weight_decay=1e-4):
    self.GAN1 = Classifier(self.model.n_latent, n_labels=self.model.n_batch, n_layers=3)
    if self.use_cuda:
        self.GAN1.cuda()
    self.optimizer_GAN = torch.optim.Adam(filter(lambda p: p.requires_grad, self.GAN1.parameters()), lr=lr,
                                          weight_decay=weight_decay)
    type(self).train(self, n_epochs=n_epochs, lr=lr)


def adversarial_wrapper(infer, warm_up=100, scale=50):
    infer.warm_up = warm_up
    infer.scale = scale
    infer.loss = types.MethodType(adversarial_loss, infer)
    infer.train = types.MethodType(adversarial_train, infer)
    return infer
