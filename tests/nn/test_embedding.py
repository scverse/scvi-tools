from __future__ import annotations

from os.path import join

import pytest
import torch

from scvi.nn import Embedding


@pytest.mark.parametrize("num_embeddings", [10])
@pytest.mark.parametrize("embedding_dim", [5])
@pytest.mark.parametrize("init", [2, [0, 1]])
@pytest.mark.parametrize("freeze_prev", [True, False])
def test_embedding_extend(
    num_embeddings: int,
    embedding_dim: int,
    init: int | list[int],
    freeze_prev: bool,
):
    embedding = Embedding(num_embeddings, embedding_dim)
    ext_embedding = Embedding.extend(embedding, init=init, freeze_prev=freeze_prev)
    n_init = len(init) if isinstance(init, list) else init

    assert ext_embedding.num_embeddings == num_embeddings + n_init
    assert ext_embedding.embedding_dim == embedding_dim
    assert ext_embedding.weight.shape == (num_embeddings + n_init, embedding_dim)
    assert torch.equal(ext_embedding.weight[:num_embeddings], embedding.weight)

    if isinstance(init, list):
        assert torch.equal(ext_embedding.weight[num_embeddings:], embedding.weight[init])

    dummy_indexes = torch.arange(num_embeddings + n_init, dtype=torch.long)
    dummy_prediction = ext_embedding(dummy_indexes)
    dummy_target = torch.randn_like(dummy_prediction)
    dummy_loss = torch.nn.functional.mse_loss(dummy_prediction, dummy_target, reduce=True)
    dummy_loss.backward()
    grad = ext_embedding.weight.grad

    if freeze_prev:
        prev_grad = grad[:num_embeddings]
        new_grad = grad[num_embeddings:]
        assert torch.equal(prev_grad, torch.zeros_like(prev_grad))
        assert not torch.equal(new_grad, torch.zeros_like(new_grad))
    else:
        assert not torch.equal(grad, torch.zeros_like(grad))


def test_embedding_extend_invalid_init(num_embeddings: int = 10, embedding_dim: int = 5):
    embedding = Embedding(num_embeddings, embedding_dim)
    with pytest.raises(ValueError):
        Embedding.extend(embedding, init=0)
    with pytest.raises(TypeError):
        Embedding.extend(embedding, init="invalid")


def test_embedding_save_load(
    save_path: str,
    num_embeddings: int = 10,
    embedding_dim: int = 5,
    init: int = 2,
    freeze_prev: bool = True,
):
    embedding = Embedding(num_embeddings, embedding_dim)
    ext_embedding = Embedding.extend(embedding, init=init, freeze_prev=freeze_prev)

    state_dict_path = join(save_path, "ext_embedding_state_dict.pt")
    torch.save(ext_embedding, state_dict_path)
    embedding = Embedding(num_embeddings, embedding_dim)
    embedding.load_state_dict(torch.load(state_dict_path, weights_only=False).state_dict())
