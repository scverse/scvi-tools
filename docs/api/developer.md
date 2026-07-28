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
   :template: class_no_inherited.rst
   :nosignatures:

   dataloaders.AnnDataLoader
   dataloaders.AnnbatchDataModule
   dataloaders.CollectionAdapter
   dataloaders.ConcatDataLoader
   dataloaders.DataSplitter
   dataloaders.SemiSupervisedDataLoader
   dataloaders.SemiSupervisedDataSplitter
   dataloaders.BatchDistributedSampler
   dataloaders.MappedCollectionDataModule
   dataloaders.TileDBDataModule

```

## Distributions

```{eval-rst}
.. currentmodule:: scvi
```

Parameterizable probability distributions.

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   distributions.Poisson
   distributions.NegativeBinomial
   distributions.NegativeBinomialMixture
   distributions.ZeroInflatedNegativeBinomial
   distributions.BetaBinomial
   distributions.Normal
   distributions.Log1pNormal
   distributions.ZeroInflatedLogNormal
   distributions.ZeroInflatedGamma

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
    model.base.SemisupervisedTrainingMixin
    model.base.PyroSviTrainMixin
    model.base.PyroSampleMixin
    model.base.PyroJitGuideWarmup
    model.base.PyroModelGuideWarmup
    model.base.DifferentialComputation
    model.base.EmbeddingMixin
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
   external.gimvi._task.GIMVITrainingPlan
   external.gimvi._task.CyclicMultiDataLoader
   external.cytovi.CytoVAE
   external.cellassign.CellAssignModule
   external.contrastivevi.ContrastiveDataSplitter
   external.contrastivevi.ContrastiveDataLoader
   external.stereoscope.RNADeconv
   external.stereoscope.SpatialDeconv
   external.scbasset.ScBassetModule
   external.contrastivevi.ContrastiveVAE
   external.velovi.VELOVAE
   external.tangram.TangramMapper
   external.mrvi.MRVAE
   external.mrvi._types.MRVIReduction
   external.methylvi.METHYLVAE
   external.methylvi.BSSeqMixin
   external.methylvi.BSSeqModuleMixin
   external.methylvi.DecoderMETHYLVI
   external.methylvi.METHYLANVAE
   external.decipher.DecipherPyroModule
   external.decipher._trainingplan.DecipherTrainingPlan
   external.resolvi.RESOLVAE
   external.totalanvi.TOTALANVAE
   external.scviva.nicheVAE
   external.scviva.NicheLossOutput
   external.scviva.differential_expression.DifferentialExpressionResults
   external.sysvi.SysVAE
   external.diagvi.DIAGVAE
   external.drvi.DRVIModule
   external.drvi.DecoderDRVI
   external.drvi.SplitFCLayers
   external.drvi.LogNegativeBinomial
   external.drvi.StackedLinearLayer
   external.JointEmbeddingVAE
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
   module.base.SupervisedModuleClass
   module.base.PyroBaseModuleClass
   module.base.EmbeddingModuleMixin
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
   nn.MultiEncoder
   nn.Decoder
   nn.MultiDecoder
   nn.DecoderSCVI
   nn.LinearDecoderSCVI
   nn.one_hot
   nn.Embedding
   nn.DecoderTOTALVI
   nn.EncoderTOTALVI

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
   train.ClassifierTrainingPlan
   train.SemiSupervisedTrainingPlan
   train.SemiSupervisedAdversarialTrainingPlan
   train.LowLevelPyroTrainingPlan
   train.PyroTrainingPlan
   train.Trainer
   train.TrainingPlan
   train.TrainRunner
   train.ScibCallback
   train.SaveCheckpoint
   train.LoudEarlyStopping
   train.KwargsConfig
   train._metrics.ElboMetric

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
   model.get_max_epochs_heuristic
   external.decipher.utils.Trajectory
```

```{eval-rst}
.. autosummary::
   :toctree: reference/
   :template: class_no_inherited.rst
   :nosignatures:

   utils.attrdict
```
