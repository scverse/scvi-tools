from __future__ import annotations

import logging
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import torch

from scvi.data import AnnDataManager
from scvi.data._download import _download
from scvi.data._preprocessing import _dna_to_code
from scvi.data.fields import CategoricalVarField, LayerField, ObsmField
from scvi.dataloaders import DataSplitter
from scvi.external.scbasset._module import REGISTRY_KEYS, ScBassetModule
from scvi.model.base import BaseModelClass
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils import dependencies, setup_anndata_dsp
from scvi.utils._docstrings import devices_dsp

if TYPE_CHECKING:
    from typing import Literal

    from anndata import AnnData

logger = logging.getLogger(__name__)


class SCBASSET(BaseModelClass):
    """``EXPERIMENTAL`` Reimplementation of scBasset :cite:p:`Yuan2022`.

    Performs representation learning of scATAC-seq data. Original implementation:
    https://github.com/calico/scBasset.

    We are working to measure the performance of this model compared to the original.

    Parameters
    ----------
    adata
        single-cell AnnData object that has been registered via
        :meth:`~scvi.external.SCBASSET.setup_anndata`.
    n_bottleneck_layer
        Size of the bottleneck layer
    l2_reg_cell_embedding
        L2 regularization for the cell embedding layer. A value, e.g. 1e-8 can be used to improve
        integration performance.
    **model_kwargs
        Keyword args for :class:`~scvi.external.scbasset.ScBassetModule`

    Examples
    --------
    >>> adata = anndata.read_h5ad(path_to_sc_anndata)
    >>> scvi.data.add_dna_sequence(adata)
    >>> adata = adata.transpose()  # regions by cells
    >>> scvi.external.SCBASSET.setup_anndata(adata, dna_code_key="dna_code")
    >>> model = scvi.external.SCBASSET(adata)
    >>> model.train()
    >>> adata.varm["X_scbasset"] = model.get_latent_representation()

    Notes
    -----
    See further usage examples in the following tutorials:

    1. :doc:`/tutorials/notebooks/atac/scbasset`
    2. :doc:`/tutorials/notebooks/atac/scbasset_batch`
    """

    MOTIF_URLS = {
        "human": (
            "https://storage.googleapis.com/scbasset_tutorial_data/Homo_sapiens_motif_fasta.tar.gz",
            "Homo_sapiens_motif_fasta",
        ),
    }
    DEFAULT_MOTIF_DIR = "./scbasset_motifs/"

    def __init__(
        self,
        adata: AnnData,
        n_bottleneck_layer: int = 32,
        l2_reg_cell_embedding: float = 0.0,
        **model_kwargs,
    ):
        super().__init__(adata)
        self.n_cells = self.summary_stats.n_vars
        self.n_regions = adata.n_obs
        self.n_batch = self.summary_stats.n_batch
        batch_ids = self.adata_manager.get_from_registry(REGISTRY_KEYS.BATCH_KEY)
        self.module = ScBassetModule(
            n_cells=self.n_cells,
            batch_ids=torch.tensor(batch_ids).long() if batch_ids.sum() > 0 else None,
            n_bottleneck_layer=n_bottleneck_layer,
            l2_reg_cell_embedding=l2_reg_cell_embedding,
            **model_kwargs,
        )
        self._model_summary_string = (
            f"ScBasset Model with params: \nn_regions: {self.n_regions}, n_batch: {self.n_batch}, "
            f"n_cells: {self.n_cells}"
        )
        self.init_params_ = self._get_init_params(locals())

    @devices_dsp.dedent
    def train(
        self,
        max_epochs: int = 1000,
        lr: float = 0.01,
        accelerator: str = "auto",
        devices: int | list[int] | str = "auto",
        train_size: float | None = None,
        validation_size: float | None = None,
        shuffle_set_split: bool = True,
        batch_size: int = 128,
        early_stopping: bool = True,
        early_stopping_monitor: str = "auroc_train",
        early_stopping_mode: Literal["min", "max"] = "max",
        early_stopping_min_delta: float = 1e-6,
        datasplitter_kwargs: dict | None = None,
        plan_kwargs: dict | None = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of epochs to train for
        lr
            Learning rate for optimization.
        %(param_accelerator)s
        %(param_devices)s
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        shuffle_set_split
            Whether to shuffle indices before splitting. If `False`, the val, train, and test set
            are split in the sequential order of the data according to `validation_size` and
            `train_size` percentages.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        early_stopping_monitor
            Metric logged during validation set epoch. The available metrics will depend on
            the training plan class used. We list the most common options here in the typing.
        early_stopping_mode
            In 'min' mode, training will stop when the quantity monitored has stopped decreasing
            and in 'max' mode it will stop when the quantity monitored has stopped increasing.
        early_stopping_min_delta
            Minimum change in the monitored quantity to qualify as an improvement,
            i.e. an absolute change of less than min_delta, will count as no improvement.
        datasplitter_kwargs
            Additional keyword arguments passed into :class:`~scvi.dataloaders.DataSplitter`.
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        custom_plan_kwargs = {
            "optimizer": "Custom",
            "optimizer_creator": lambda p: torch.optim.Adam(p, lr=lr, betas=(0.95, 0.9995)),
        }
        if plan_kwargs is not None:
            custom_plan_kwargs.update(plan_kwargs)

        datasplitter_kwargs = datasplitter_kwargs or {}

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            shuffle_set_split=shuffle_set_split,
            batch_size=batch_size,
            # We don't want to dataload the batch ids into the module
            data_and_attributes={
                REGISTRY_KEYS.X_KEY: np.float32,
                REGISTRY_KEYS.DNA_CODE_KEY: np.int64,
            },
            **datasplitter_kwargs,
        )
        training_plan = TrainingPlan(self.module, **custom_plan_kwargs)

        es = {
            "early_stopping": early_stopping,
            "early_stopping_monitor": early_stopping_monitor,
            "early_stopping_mode": early_stopping_mode,
            "early_stopping_min_delta": early_stopping_min_delta,
        }
        for k, v in es.items():
            trainer_kwargs[k] = v if k not in trainer_kwargs.keys() else trainer_kwargs[k]
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            accelerator=accelerator,
            devices=devices,
            **trainer_kwargs,
        )
        return runner()

    @torch.inference_mode()
    def get_latent_representation(self) -> np.ndarray:
        """Returns the latent representation of the cells.

        Returns
        -------
        latent representation (n_cells, n_latent)
        """
        return self.module.cell_embedding.cpu().numpy().T

    @torch.inference_mode()
    def get_cell_bias(self) -> np.ndarray:
        """Returns the cell-specific bias term.

        Returns
        -------
        bias (n_cells,)
        """
        return self.module.cell_bias.cpu().numpy()

    @classmethod
    def _download_motifs(cls, genome: str, motif_dir: str) -> None:
        """Download a set of motifs injected into peaks"""
        logger.info(f"Downloading motif set to: {motif_dir}")
        # download the motif set
        url_name = cls.MOTIF_URLS.get(genome, None)  # (url, dir_name)
        if url_name is None:
            raise ValueError(f"{genome} is not a supported motif set.")
        _download(url_name[0], save_path=motif_dir, filename=f"{genome}_motifs.tar.gz")
        # untar it
        import tarfile

        def rename_members(tarball):
            """Rename files in the tarball to remove the top level folder"""
            for member in tarball.getmembers():
                if member.path.startswith(url_name[1]):
                    member.path = member.path.replace(url_name[1] + "/", "")
                    yield member

        # importing the "tarfile" module
        tarball = tarfile.open(Path(motif_dir, f"{genome}_motifs.tar.gz"))
        tarball.extractall(path=motif_dir, members=rename_members(tarball))
        tarball.close()

        # `motif_dir` now has `shuffled_peaks_motifs` as a subdir and
        # `shuffled_peaks.fasta` as a root level file.
        logger.info("Download and extraction complete.")
        return

    @dependencies("Bio")
    def _get_motif_library(
        self, tf: str, genome: str = "human", motif_dir: str | None = None
    ) -> tuple[list[str], list[str]]:
        """Load sequences with a TF motif injected from a pre-computed library

        Parameters
        ----------
        tf
            name of the transcription factor motif to load. Must be present in a
            pre-computed library.
        genome
            species name for the motif injection procedure. Currently, only "human"
            is supported.
        motif_dir
            path for the motif library. Will download if not already present.

        Returns
        -------
        motif_seqs
            list of sequences with an injected motif.
        bg_seqs
            dinucleotide shuffled background sequences.
        """
        from Bio import SeqIO

        if motif_dir is None:
            motif_dir = self.DEFAULT_MOTIF_DIR

        # ensure input is a `Path` object
        motif_dir = Path(motif_dir)
        if not Path(motif_dir, "shuffled_peaks.fasta").exists():
            motif_dir.mkdir(exist_ok=True, parents=True)
            self._download_motifs(genome=genome, motif_dir=motif_dir)

        fasta_files = motif_dir.glob("shuffled_peaks_motifs/*.fasta")
        tf_names = [f.stem for f in fasta_files]
        if tf not in tf_names:
            msg = f"{tf} is not found as a motif in the library."
            raise ValueError(msg)

        # load the motifs
        tf_motif_path = Path(motif_dir, "shuffled_peaks_motifs", f"{tf}.fasta")
        motif_seqs = list(SeqIO.parse(tf_motif_path, "fasta"))
        motif_seqs = [str(i.seq) for i in motif_seqs]

        bg_seqs_path = Path(motif_dir, "shuffled_peaks.fasta")
        bg_seqs = list(SeqIO.parse(bg_seqs_path, "fasta"))
        bg_seqs = [str(i.seq) for i in bg_seqs]
        return motif_seqs, bg_seqs

    @torch.inference_mode()
    def get_tf_activity(
        self,
        tf: str,
        genome: str = "human",
        motif_dir: str | None = None,
        lib_size_norm: bool | None = True,
        batch_size: int = 256,
    ) -> np.ndarray:
        """Infer transcription factor activity using a motif injection procedure.

        Parameters
        ----------
        tf
            transcription factor name. must be provided in the relevant motif repository.
        genome
            species name for the motif injection procedure. Currently, only "human"
            is supported.
        motif_dir
            path for the motif library. Will download if not already present.
        lib_size_norm
            normalize accessibility scores for library size by *substracting* the
            cell bias term from each accessibility score prior to comparing motif
            scores to background scores.
        batch_size
            minibatch size for TF activity inference.

        Returns
        -------
        tf_score
            [cells,] TF activity scores.

        Notes
        -----
        scBasset infers TF activities by injecting known TF motifs into a
        shuffled dinucleotide sequence and computing the change in accessibility
        predicted between the injected motif and a randomized background
        sequence. See :cite:p:`Yuan2022` for details. We modeled this function
        off the original implementation in `scbasset`.

        https://github.com/calico/scBasset/blob/9a3294c240d849cdac62682e324bc5f4836bb744/scbasset/utils.py#L453
        """
        # check for a library of FASTA sequences corresponding to motifs and
        # download if none is found
        # `motif_seqs` is a List of str sequences where each char is in "ACTGN".
        # `bg_seqs` is the same, but background sequences rather than motif injected
        motif_seqs, bg_seqs = self._get_motif_library(tf=tf, genome=genome, motif_dir=motif_dir)

        # SCBASSET.module.inference(...) takes `dna_code: torch.Tensor` as input
        # where `dna_code` is [batch_size, seq_length] and each value is [0,1,2,3]
        # where [0: A, 1: C, 2: G, 3: T].
        motif_codes = pd.DataFrame([list(s) for s in motif_seqs]).applymap(_dna_to_code)
        bg_codes = pd.DataFrame([list(s) for s in bg_seqs]).applymap(_dna_to_code)
        # [batch_size, seq_length]
        motif_codes = torch.from_numpy(np.array(motif_codes)).long()
        bg_codes = torch.from_numpy(np.array(bg_codes)).long()

        # NOTE: SCBASSET uses a fixed size of 1344 bp. If motifs from a different source
        # than the above are used, we may need to truncate to match the model size.
        # We should be cautious about doing this, so we throw a warning to the user.
        model_input_size = self.adata_manager.get_from_registry(REGISTRY_KEYS.DNA_CODE_KEY).shape[
            1
        ]
        n_diff = motif_codes.shape[1] - model_input_size
        if n_diff > 0:
            n_cut = n_diff // 2
            logger.warning(
                f"Motif size {motif_codes.shape[1]} != model input size {model_input_size}."
                f" Trimming {n_cut} from each side."
            )
            motif_codes = motif_codes[:, n_cut:-n_cut]
            bg_codes = bg_codes[:, n_cut:-n_cut]
        if n_diff < 0:
            msg = f"Motif sizes {motif_codes.shape[1]} < model size {model_input_size}"
            raise ValueError(msg)

        motif_accessibility = self.module._get_accessibility(
            dna_codes=motif_codes,
            batch_size=batch_size,
        )
        bg_accessibility = self.module._get_accessibility(
            dna_codes=bg_codes,
            batch_size=batch_size,
        )
        # move to CPU
        motif_accessibility = motif_accessibility.detach().cpu()
        bg_accessibility = bg_accessibility.detach().cpu()
        if lib_size_norm:
            # substract the cell bias term so that scores are agnostic to the
            # library size of each observation
            bias = self.module.cell_bias.detach().cpu()
            motif_accessibility = motif_accessibility - bias
            bg_accessibility = bg_accessibility - bias

        # compute the difference in activity between the motif and background
        # sequences
        # after means, arr is activity by cell, shape [cells,]
        motif_activity = motif_accessibility.mean(0) - bg_accessibility.mean(0)
        motif_activity = motif_activity.numpy()
        # z-score the activity
        tf_score = (motif_activity - motif_activity.mean()) / motif_activity.std()
        return tf_score

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        dna_code_key: str,
        layer: str | None = None,
        batch_key: str | None = None,
        **kwargs,
    ):
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        dna_code_key
            Key in `adata.obsm` with dna sequences encoded as integer code.
        %(param_layer)s
        batch_key
            key in `adata.var` for batch information. Categories will automatically be converted
            into integer categories and saved to `adata.var['_scvi_batch']`. If `None`, assigns the
            same batch to all the data.

        Notes
        -----
        The adata object should be in the regions by cells format. This is due to scBasset
        considering regions as observations and cells as variables. This can be simply achieved
        by transposing the data, `bdata = adata.transpose()`.
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, layer, is_count_data=True),
            ObsmField(REGISTRY_KEYS.DNA_CODE_KEY, dna_code_key, is_count_data=True),
            CategoricalVarField(REGISTRY_KEYS.BATCH_KEY, batch_key),
        ]
        adata_manager = AnnDataManager(fields=anndata_fields, setup_method_args=setup_method_args)
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)
