import torch

import scvi
from scvi import REGISTRY_KEYS


def test_contrastive_dataloader(
    mock_contrastive_adata_manager,
    mock_background_indices,
    mock_target_indices,
):
    batch_size = 32
    adata = mock_contrastive_adata_manager.adata

    expected_background_indices = mock_background_indices[:batch_size]
    expected_background_input = torch.Tensor(
        adata.layers["raw_counts"][expected_background_indices, :]
    )
    expected_background_labels = torch.LongTensor(
        adata.obs["_scvi_labels"].iloc[expected_background_indices]
    ).unsqueeze(1)
    expected_background_batch = torch.LongTensor(
        adata.obs["_scvi_batch"].iloc[expected_background_indices]
    ).unsqueeze(1)

    expected_target_indices = mock_target_indices[:batch_size]
    expected_target_input = torch.Tensor(adata.layers["raw_counts"][expected_target_indices, :])
    expected_target_labels = torch.LongTensor(
        adata.obs["_scvi_labels"].iloc[expected_target_indices]
    ).unsqueeze(1)
    expected_target_batch = torch.LongTensor(
        adata.obs["_scvi_batch"].iloc[expected_target_indices]
    ).unsqueeze(1)

    dataloader = scvi.external.contrastivevi.ContrastiveDataLoader(
        adata_manager=mock_contrastive_adata_manager,
        background_indices=mock_background_indices,
        target_indices=mock_target_indices,
        batch_size=batch_size,
    )
    batch = next(batch for batch in dataloader)

    assert isinstance(batch, dict)
    assert len(batch.keys()) == 2
    assert "background" in batch.keys()
    assert "target" in batch.keys()

    assert torch.equal(batch["background"][REGISTRY_KEYS.X_KEY], expected_background_input)
    assert torch.equal(batch["background"][REGISTRY_KEYS.LABELS_KEY], expected_background_labels)
    assert torch.equal(batch["background"][REGISTRY_KEYS.BATCH_KEY], expected_background_batch)

    assert torch.equal(batch["target"][REGISTRY_KEYS.X_KEY], expected_target_input)
    assert torch.equal(batch["target"][REGISTRY_KEYS.LABELS_KEY], expected_target_labels)
    assert torch.equal(batch["target"][REGISTRY_KEYS.BATCH_KEY], expected_target_batch)
