from typing import NamedTuple


class _SCVI_HUB_NT(NamedTuple):
    HF_LIBRARY_NAME: str = "scvi-tools"
    MAX_HF_UPLOAD_SIZE: int = 5e9  # 5GB

    # file names
    METADATA_FILE_NAME: str = "_scvi_required_metadata.json"
    MODEL_CARD_FILE_NAME: str = "README.md"

    # model card defaults
    DEFAULT_MISSING_FIELD: str = "To be added..."
    DEFAULT_NA_FIELD: str = "N/A"
    DEFAULT_PARENT_MODULE: str = "scvi.model"

    # model card tags
    MODEL_CLS_NAME_TAG: str = "model_cls_name:{}"
    SCVI_VERSION_TAG: str = "scvi_version:{}"
    ANNDATA_VERSION_TAG: str = "anndata_version:{}"
    MODALITY_TAG: str = "modality:{}"
    TISSUE_TAG: str = "tissue:{}"
    ANNOTATED_TAG: str = "annotated:{}"


_SCVI_HUB = _SCVI_HUB_NT()
