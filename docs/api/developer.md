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
   data.AnnDataManagerValidationCheck
   data.fields.BaseAnnDataField
   data.fields.LayerField
   data.fields.CategoricalObsField
   data.fields.CategoricalVarField
   data.fields.NumericalJointObsField
   data.fields.NumericalJointVarField
   data.fields.CategoricalJointObsField
   data.fields.CategoricalJointVarField
   data.fields.ObsmField
   data.fields.VarmField
   data.fields.ProteinObsmField
   data.fields.StringUnsField
   data.fields.LabelsWithUnlabeledObsField
   data.fields.BaseMuDataWrapperClass
   data.fields.MuDataWrapper
   data.fields.MuDataLayerField
   data.fields.MuDataProteinLayerField
   data.fields.MuDataNumericalObsField
   data.fields.MuDataNumericalVarField
   data.fields.MuDataCategoricalObsField
   data.fields.MuDataCategoricalVarField
   data.fields.MuDataObsmField
   data.fields.MuDataVarmField
   data.fields.MuDataNumericalJointObsField
   data.fields.MuDataNumericalJointVarField
   data.fields.MuDataCategoricalJointObsField
   data.fields.MuDataCategoricalJointVarField
   data.AnnTorchDataset

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
    model.base.BaseMinifiedModeModelClass
    model.base.VAEMixin
    model.base.RNASeqMixin
    model.base.ArchesMixin
    model.base.UnsupervisedTrainingMixin
    model.base.PyroSviTrainMixin
    model.base.PyroSampleMixin
    model.base.PyroJitGuideWarmup
    model.base.PyroModelGuideWarmup
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
   module.JaxVAE

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
   external.tangram.TangramMapper
   external.scbasset.ScBassetModule

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
   module.base.BaseMinifiedModeModuleClass
   module.base.PyroBaseModuleClass
   module.base.JaxBaseModuleClass
   module.base.LossOutput
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
   train.SemiSupervisedTrainingPlan
   train.LowLevelPyroTrainingPlan
   train.PyroTrainingPlan
   train.JaxTrainingPlan
   train.Trainer
   train.TrainingPlan
   train.TrainRunner
   train.SaveBestState
   train.LoudEarlyStopping

```

## Model hyperparameter tuning

`scvi-tools` supports automatic model hyperparameter tuning using [Ray Tune]. These
classes allow for new model classes to be easily integrated with the module.

```{eval-rst}
.. currentmodule:: scvi
```

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :nosignatures:

   autotune.TunerManager
   autotune.Tunable
   autotune.TunableMixin
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

[ray tune]: https://docs.ray.io/en/latest/tune/index.html
