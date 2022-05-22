from typing import Iterable, Optional, Sequence

import numpy as np
import torch
from torch.distributions import Categorical, Normal  #ok
from torch.distributions import kl_divergence as kl
from torch.nn import functional as F

from scvi import REGISTRY_KEYS
from scvi._compat import Literal
from scvi.module.base import LossRecorder, auto_move_data
from scvi.nn import Decoder, Encoder

from ._classifier import Classifier   #Basic fully-connected NN classifier.
from ._utils import broadcast_labels
from ._vae import VAE


class SCANVAE(VAE): #inherits from VAE class (for instance inherits z_encoder)
    """
    Single-cell annotation using variational inference.

    This is an implementation of the scANVI model described in [Xu21]_,
    inspired from M1 + M2 model, as described in (https://arxiv.org/pdf/1406.5298.pdf).

    Parameters
    ----------
    n_version
    n_input
        Number of input genes
    n_batch
        Number of batches
    n_labels
        Number of labels
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    n_continuous_cov
        Number of continuous covariates  
    n_cats_per_cov
        Number of categories for each extra categorical covariate
    dropout_rate
        Dropout rate for neural networks
    dispersion
        One of the following

        * ``'gene'`` - dispersion parameter of NB is constant per gene across cells
        * ``'gene-batch'`` - dispersion can differ between different batches
        * ``'gene-label'`` - dispersion can differ between different labels
        * ``'gene-cell'`` - dispersion can differ for every gene in every cell
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    gene_likelihood
        One of

        * ``'nb'`` - Negative binomial distribution
        * ``'zinb'`` - Zero-inflated negative binomial distribution
    y_prior
        If None, initialized to uniform probability over cell types  OK
    labels_groups
        Label group designations                                                    ?? --> hierarchie entre labels
    use_labels_groups
        Whether to use the label groups
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    **vae_kwargs
        Keyword args for :class:`~scvi.module.VAE`
    """


    #--------------------------------INIT-----------------------------------------------------------------------------------------------------------

    def __init__(
        self,
        n_input: int,
        n_batch: int = 0,
        n_labels: int = 0,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        n_continuous_cov: int = 0, #in the following, we assume only one categorical covariate with categories, which represents the common case of having multiple batches of data.

        n_cats_per_cov: Optional[Iterable[int]] = None,
        dropout_rate: float = 0.1,
        dispersion: str = "gene",
        log_variational: bool = True,
        gene_likelihood: str = "zinb",
        y_prior=None,
        labels_groups: Sequence[int] = None,                #??
        use_labels_groups: bool = False,
        classifier_parameters: dict = dict(),
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "none",
        n_version = 0,  # 0 denotes the old one without the fix
        **vae_kwargs
    ):
        super().__init__(
            n_input,
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            n_continuous_cov=n_continuous_cov,
            n_cats_per_cov=n_cats_per_cov,
            dropout_rate=dropout_rate,
            n_batch=n_batch,
            dispersion=dispersion,
            log_variational=log_variational,
            gene_likelihood=gene_likelihood,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **vae_kwargs
        )

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"

        self.n_version = n_version
        self.n_labels = n_labels


        # Classifier takes n_latent as input
        cls_parameters = {
            "n_layers": n_layers,
            "n_hidden": n_hidden,
            "dropout_rate": dropout_rate,
        }
        cls_parameters.update(classifier_parameters)
        self.classifier = Classifier(                           #PROBABILISTIC CELL-TYPE ANNOTATION? n_hidden kept as default, classifies between n_labels
            n_latent,  #Number of input dimensions
            n_labels=n_labels,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            **cls_parameters
        )

        self.encoder_z2_z1 = Encoder(                #q(z2|z1,....)  ???
            n_latent,
            n_latent,
            n_cat_list=[self.n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
        )
        self.decoder_z1_z2 = Decoder(               #  p(z1|z2,....)  ????
            n_latent,
            n_latent,
            n_cat_list=[self.n_labels],
            n_layers=n_layers,
            n_hidden=n_hidden,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
        )

        self.y_prior = torch.nn.Parameter(    #uniform probabilities for categorical distribution on the cell type  HERE Y=C
            y_prior
            if y_prior is not None
            else (1 / n_labels) * torch.ones(1, n_labels),
            requires_grad=False,
        )
        self.use_labels_groups = use_labels_groups
        self.labels_groups = (
            np.array(labels_groups) if labels_groups is not None else None
        )
        if self.use_labels_groups:           
            if labels_groups is None:
                raise ValueError("Specify label groups")
            unique_groups = np.unique(self.labels_groups)
            self.n_groups = len(unique_groups)
            if not (unique_groups == np.arange(self.n_groups)).all():
                raise ValueError()

            self.classifier_groups = Classifier(
                n_latent, n_hidden, self.n_groups, n_layers, dropout_rate
            )
            self.groups_index = torch.nn.ParameterList(
                [
                    torch.nn.Parameter(
                        torch.tensor(
                            (self.labels_groups == i).astype(np.uint8),
                            dtype=torch.uint8,
                        ),
                        requires_grad=False,
                    )
                    for i in range(self.n_groups)
                ]
            )


    #---------------------------------------METHODS----------------------------------------------------------------------------------------------------------------------------

    @auto_move_data
    def classify(self, x, batch_index=None, cont_covs=None, cat_covs=None):
        if self.log_variational:   #for numerical stability
            x = torch.log(1 + x)

        if cont_covs is not None and self.encode_covariates:
            encoder_input = torch.cat((x, cont_covs), dim=-1)
        else:
            encoder_input = x
        if cat_covs is not None and self.encode_covariates:
            categorical_input = torch.split(cat_covs, 1, dim=1)
        else:
            categorical_input = tuple()

        qz_m, _, z = self.z_encoder(encoder_input, batch_index, *categorical_input)    #q(z1|x)  without the var   qz_v 
        # We classify using the inferred mean parameter of z_1 in the latent space
        z = qz_m    
        if self.use_labels_groups:
            w_g = self.classifier_groups(z)
            unw_y = self.classifier(z)
            w_y = torch.zeros_like(unw_y)
            for i, group_index in enumerate(self.groups_index):
                unw_y_g = unw_y[:, group_index]
                w_y[:, group_index] = unw_y_g / (
                    unw_y_g.sum(dim=-1, keepdim=True) + 1e-8
                )
                w_y[:, group_index] *= w_g[:, [i]]
        else:
            w_y = self.classifier(z)
        return w_y

    @auto_move_data
    def classification_loss(self, labelled_dataset):       #add a classifiaction loss ON THE LABELLED ATA, following Kingma et al 
        x = labelled_dataset[REGISTRY_KEYS.X_KEY]
        y = labelled_dataset[REGISTRY_KEYS.LABELS_KEY]
        batch_idx = labelled_dataset[REGISTRY_KEYS.BATCH_KEY]
        cont_key = REGISTRY_KEYS.CONT_COVS_KEY
        cont_covs = (
            labelled_dataset[cont_key] if cont_key in labelled_dataset.keys() else None
        )

        cat_key = REGISTRY_KEYS.CAT_COVS_KEY
        cat_covs = (
            labelled_dataset[cat_key] if cat_key in labelled_dataset.keys() else None
        )
        classification_loss = F.cross_entropy(
            self.classify(
                x, batch_index=batch_idx, cat_covs=cat_covs, cont_covs=cont_covs
            ),
            y.view(-1).long(),  
        )
        return classification_loss

    def loss(
        self,
        tensors,
        inference_outputs,
        generative_ouputs,
        feed_labels=False,  #? ---> 2 dataloaders, for annotated and un annotated, don't feed labels for un annotated
        kl_weight=1,
        labelled_tensors=None,  #??  -->scvanvi.py
        classification_ratio=None,
    ):
        px_r = generative_ouputs["px_r"]
        px_rate = generative_ouputs["px_rate"]
        px_dropout = generative_ouputs["px_dropout"]
        qz1_m = inference_outputs["qz_m"]
        qz1_v = inference_outputs["qz_v"]
        z1 = inference_outputs["z"]
        x = tensors[REGISTRY_KEYS.X_KEY]
        batch_index = tensors[REGISTRY_KEYS.BATCH_KEY]

        if feed_labels:
            y = tensors[REGISTRY_KEYS.LABELS_KEY]
        else:
            y = None


        is_labelled = False if y is None else True  #important for ELBO 

        # Enumerate choices of label
        ys, z1s = broadcast_labels(y, z1, n_broadcast=self.n_labels) #one-hot encoding of the labels
        #if z1 is of size (batch_size,latent), z1_s is of size (n_labels*batch_size,latent) 
        qz2_m, qz2_v, z2 = self.encoder_z2_z1(z1s, ys)  #q(z2|z1,..)
        pz1_m, pz1_v = self.decoder_z1_z2(z2, ys)  #p(z1|z2,..)

        reconst_loss = self.get_reconstruction_loss(x, px_rate, px_r, px_dropout)  #expectation of log(p),as in scvi

        # KL Divergence
        mean = torch.zeros_like(qz2_m)
        scale = torch.ones_like(qz2_v)

        kl_divergence_z2 = kl(
            Normal(qz2_m, torch.sqrt(qz2_v)), Normal(mean, scale)  #q(z2|z1,..)||p(z2)
        ).sum(dim=1)

        loss_z1_unweight = -Normal(pz1_m, torch.sqrt(pz1_v)).log_prob(z1s).sum(dim=-1)
        #Sum of the log of the  Normal probability density evaluated at value z1s. The sum is over the 10-dim latent space.


        loss_z1_weight = Normal(qz1_m, torch.sqrt(qz1_v)).log_prob(z1).sum(dim=-1)  #????ok
        
        if not self.use_observed_lib_size:   #comme dans le user guide, si l_n is latent!
            ql_m = inference_outputs["ql_m"]
            ql_v = inference_outputs["ql_v"]
            (
                local_library_log_means,
                local_library_log_vars,
            ) = self._compute_local_library_params(batch_index)

            kl_divergence_l = kl(
                Normal(ql_m, torch.sqrt(ql_v)),
                Normal(local_library_log_means, torch.sqrt(local_library_log_vars)),   #ok
            ).sum(dim=1)
        else:
            kl_divergence_l = 0.0  #indeed si tu observes l il ne sera plus dans la var dist

        print('The version is: ', self.n_version)

        #if is_labelled:
        if labelled_tensors is not None:
            print("--------------------labelled_tensors is not None-------------------------")
            if self.n_version == 1:
                print("Adding KLs to the loss...")
                loss = reconst_loss + loss_z1_weight + loss_z1_unweight + kl_divergence_z2 + kl_divergence_l  # add kl terms here
            else:
                print("The loss is unchanged...")
                loss = reconst_loss + loss_z1_weight + loss_z1_unweight

            kl_locals = {
                "kl_divergence_z2": kl_divergence_z2,  #in scvi, this is added to the loss?
                "kl_divergence_l": kl_divergence_l,
            }
            #if labelled_tensors is not None:
                #print("And labelled_tensors is not None")
            classifier_loss = self.classification_loss(labelled_tensors)
            loss += classifier_loss * classification_ratio
            return LossRecorder(
                loss,
                reconst_loss,
                kl_locals,
                classification_loss=classifier_loss,
                n_labelled_tensors=labelled_tensors[REGISTRY_KEYS.X_KEY].shape[0],
            )
            #return LossRecorder(
            #    loss,
            #    reconst_loss,
            #    kl_locals,
            #    kl_global=torch.tensor(0.0),
            #)

        # the ELBO in the case where C=Y is not observed
        probs = self.classifier(z1)  #outputs a vector of size n_labels suming to 1
        reconst_loss += loss_z1_weight + (
            (loss_z1_unweight).view(self.n_labels, -1).t() * probs
        ).sum(dim=1)  #why loss_z1_weight is not in the sum? 

        kl_divergence = (kl_divergence_z2.view(self.n_labels, -1).t() * probs).sum(
            dim=1
        )
        kl_divergence += kl(
            Categorical(probs=probs),
            Categorical(probs=self.y_prior.repeat(probs.size(0), 1)),
        )
        kl_divergence += kl_divergence_l

        loss = torch.mean(reconst_loss + kl_divergence * kl_weight)   #annealing to avoid posterior collapse!!!

       # if labelled_tensors is not None:
       #     print("is_labelled=False and labelled_tensors is not None")
       #     classifier_loss = self.classification_loss(labelled_tensors)
       #     loss += classifier_loss * classification_ratio
       #     return LossRecorder(
       #         loss,
       #         reconst_loss,
       #         kl_divergence,
       #         classification_loss=classifier_loss,
       #     )
       return LossRecorder(loss, reconst_loss, kl_divergence)
