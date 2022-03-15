# Developer

Import scvi-tools as:

```
import scvi
```

```{eval-rst}
.. currentmodule:: scvi
```

## Data Registration

```{eval-rst}
.. currentmodule:: scvi
```

AnnDataFields delineate how scvi-tools refers to fields in AnnData objects. The AnnDataManager provides an interface
for operating over a collection of AnnDataFields and an AnnData object.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   data.AnnDataManager
   data.fields.BaseAnnDataField
   data.fields.LayerField
   data.fields.CategoricalObsField
   data.fields.NumericalJointObsField
   data.fields.CategoricalJointObsField
   data.fields.ObsmField
   data.fields.ProteinObsmField
   data.fields.LabelsWithUnlabeledObsField

```

## Data Loaders

```{eval-rst}
.. currentmodule:: scvi
```

DataLoaders for loading tensors from AnnData objects. DataSplitters for splitting data into train/test/val.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   dataloaders.AnnDataLoader
   dataloaders.AnnTorchDataset
   dataloaders.ConcatDataLoader
   dataloaders.DataSplitter
   dataloaders.SemiSupervisedDataLoader
   dataloaders.SemiSupervisedDataSplitter

```

## Distributions

```{eval-rst}
.. currentmodule:: scvi
```

Parameterizable probability distributions.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   distributions.NegativeBinomial
   distributions.NegativeBinomialMixture
   distributions.ZeroInflatedNegativeBinomial
   distributions.JaxNegativeBinomialMeanDisp

```

## Model (Base)

```{eval-rst}
.. currentmodule:: scvi
```

These classes should be used to construct user-facing model classes.

```{eval-rst}
.. autosummary::
    :toctree: reference/
    :nosignatures:

    model.base.BaseModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin
    model.base.UnsupervisedTrainingMixin
    model.base.PyroSviTrainMixin
    model.base.PyroSampleMixin
    model.base.PyroJitGuideWarmup
    model.base.DifferentialComputation
```

## Module

```{eval-rst}
.. currentmodule:: scvi
```

Existing module classes with respective generative and inference procedures.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   module.AutoZIVAE
   module.Classifier
   module.LDVAE
   module.MRDeconv
   module.PEAKVAE
   module.MULTIVAE
   module.SCANVAE
   module.TOTALVAE
   module.VAE
   module.VAEC
   module.AmortizedLDAPyroModule

```

## External module

```{eval-rst}
.. currentmodule:: scvi
```

Module classes in the external API with respective generative and inference procedures.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   external.gimvi.JVAE
   external.cellassign.CellAssignModule
   external.stereoscope.RNADeconv
   external.stereoscope.SpatialDeconv

```

## Module (Base)

```{eval-rst}
.. currentmodule:: scvi
```

These classes should be used to construct module classes that define generative models and inference schemes.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   module.base.BaseModuleClass
   module.base.PyroBaseModuleClass
   module.base.LossRecorder
   module.base.auto_move_data

```

## Neural networks

```{eval-rst}
.. currentmodule:: scvi
```

Basic neural network building blocks.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   nn.FCLayers
   nn.Encoder
   nn.Decoder
   nn.one_hot

```

## Train

```{eval-rst}
.. currentmodule:: scvi

```

TrainingPlans define train/test/val optimization steps for modules.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   train.AdversarialTrainingPlan
   train.PyroTrainingPlan
   train.SemiSupervisedTrainingPlan
   train.Trainer
   train.TrainingPlan
   train.TrainRunner
   train.SaveBestState
   train.LoudEarlyStopping

```

## Utilities

```{eval-rst}
.. currentmodule:: scvi
```

Utility functions used by scvi-tools.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   utils.track
   utils.setup_anndata_dsp
   utils.attrdict
```
