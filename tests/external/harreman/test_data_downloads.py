from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

import scvi.external.harreman.datasets.datasets as datasets_mod
import scvi.external.harreman.hotspot.local_autocorrelation as autocorr_mod
import scvi.external.harreman.preprocessing.database as database_mod

BASE_URL = "https://exampledata.scverse.org/scvi-tools/harreman"

DATASET_HASHES = {
    "Liu_et_al_human_lung.h5ad": (
        "sha256:e75a13e6931194968c5ec83a986acf63addcefde8bbfc13039b76e5112d52d5e"
    ),
    "Liu_et_al_human_lung_Puck_200727_08.h5ad": (
        "sha256:a221f88cc5473f371768ac3c41d80b602f865e1b2644183829a8a595570a0d51"
    ),
    "Liu_et_al_human_lung_Puck_200727_09.h5ad": (
        "sha256:c6b8d27c23c21373542c7501fa84d4543eb2abd0fab9101ade69264059c3902a"
    ),
    "Liu_et_al_human_lung_Puck_200727_10.h5ad": (
        "sha256:54df57d215ad7d0c465e50f963b467b9b20c2970a94394df49489f9a86628dcc"
    ),
    "Liu_et_al_human_lung_Puck_220408_13.h5ad": (
        "sha256:78bfe2e057ea8d35c03eeb81731979abc07c9c00f2707c7b18e388c71541fa82"
    ),
    "Liu_et_al_human_lung_Puck_220408_14.h5ad": (
        "sha256:6e6c4df98b52cbbe386124626c45b5134a5af63be01607339263f882dd78f5a8"
    ),
    "Liu_et_al_human_lung_Puck_220408_15.h5ad": (
        "sha256:c3d0101de3911c009ad0897e0dbd7ceb430b58ecc44218f12fb7892f6c81f100"
    ),
    "Liu_et_al_human_lung_Puck_220408_20.h5ad": (
        "sha256:2cfcbd0331d4730f34640284fb8324654fd5856a48f6488521c7e71a06d9973b"
    ),
    "Parigi_et_al_mouse_colon_d0.h5ad": (
        "sha256:32467aca0f5dd4781e3c231b1c056d0156230b5df946b5f1d6028caabff487fc"
    ),
    "Parigi_et_al_mouse_colon_d14.h5ad": (
        "sha256:86b0bc9613c673ec5d153ed78d7150443691d02fd6a4f00e50314071e093d66c"
    ),
    "Parigi_et_al_mouse_colon_unrolled.h5ad": (
        "sha256:de312840a7fcd3252cc558906e505c4dcfbd7159fa658669c465c2efc50b10db"
    ),
}

HARREMAN_DB_HASHES = {
    "HarremanDB_human.csv": (
        "sha256:34d389de1d3a039b2f481444ba9f868da249c771e6d46ecdc89712be95b19258"
    ),
    "HarremanDB_human_extracellular.csv": (
        "sha256:9848602df48fb70c59060ac378016d4eb58c826ddd754af427379d43ef1be789"
    ),
    "HarremanDB_mouse.csv": (
        "sha256:d05c1eeda6975a9a25c5b007ddb696ff41f43bca525103a7164196bd54b92f24"
    ),
    "HarremanDB_mouse_extracellular.csv": (
        "sha256:7bfbb2b4c10d5e9ba8b5f2cb4d3c86cc2252c657887820a611e0e02b424850a3"
    ),
    "Heterodimer_info_human.csv": (
        "sha256:e2e5f53a0fcdeebedd012c9c9924a8a4ff451e3ba0658d37f56f11450651b063"
    ),
    "Heterodimer_info_mouse.csv": (
        "sha256:1476262cd088c64f3e7852358082f776867c8a9231f0b01335516424d07253c9"
    ),
}

CELLCHAT_HASHES = {
    "complex_input_CellChatDB_v2_human.csv": (
        "sha256:3d683825466f69a8f4f3bfb323028168e5fd2bc99fd4641b1ad10271a27aad40"
    ),
    "complex_input_CellChatDB_v2_mouse.csv": (
        "sha256:90e285fd737ed5919e8f7debbbc1b7270350997c16026da38c1360d229a1ed2e"
    ),
    "interaction_input_CellChatDB_v2_human.csv": (
        "sha256:4e04eeb9afc062c1742df26cffde68005ffa5931b86187f7a50ef583b3ea7cb5"
    ),
    "interaction_input_CellChatDB_v2_mouse.csv": (
        "sha256:48b3b493c624b267b88a8dc709dbac81c3e6c601a5d0346f3073ce50a0abb1bc"
    ),
}

METABOLIC_GENE_HASHES = {
    "human_metabolic_genes.csv": (
        "sha256:e9b726a2d06f63e822c3ed323d5a67628e12c7921eb46599a9a3a718849eb807"
    ),
    "mouse_metabolic_genes.csv": (
        "sha256:cc0093c42fab75a3a3310cd12f63ede6a5a0d400a5ca679114b11cbba8b1f1aa"
    ),
}


def _patch_retrieve(monkeypatch, module):
    calls = []

    def retrieve(**kwargs):
        calls.append(kwargs)
        return f"/tmp/harreman-cache/{kwargs['fname']}"

    monkeypatch.setattr(module.pooch, "retrieve", retrieve)
    return calls


@pytest.mark.parametrize(
    ("loader", "sample", "fname"),
    [
        (
            datasets_mod.load_visium_mouse_colon_dataset,
            None,
            "Parigi_et_al_mouse_colon_unrolled.h5ad",
        ),
        (
            datasets_mod.load_visium_mouse_colon_dataset,
            "d0",
            "Parigi_et_al_mouse_colon_d0.h5ad",
        ),
        (
            datasets_mod.load_visium_mouse_colon_dataset,
            "d14",
            "Parigi_et_al_mouse_colon_d14.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            None,
            "Liu_et_al_human_lung.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_200727_08",
            "Liu_et_al_human_lung_Puck_200727_08.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_200727_09",
            "Liu_et_al_human_lung_Puck_200727_09.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_200727_10",
            "Liu_et_al_human_lung_Puck_200727_10.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_220408_13",
            "Liu_et_al_human_lung_Puck_220408_13.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_220408_14",
            "Liu_et_al_human_lung_Puck_220408_14.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_220408_15",
            "Liu_et_al_human_lung_Puck_220408_15.h5ad",
        ),
        (
            datasets_mod.load_slide_seq_human_lung_dataset,
            "Puck_220408_20",
            "Liu_et_al_human_lung_Puck_220408_20.h5ad",
        ),
    ],
)
def test_dataset_loaders_fetch_from_exampledata_with_hash(monkeypatch, loader, sample, fname):
    calls = _patch_retrieve(monkeypatch, datasets_mod)
    read_paths = []
    sentinel = object()

    def read_h5ad(path):
        read_paths.append(path)
        return sentinel

    monkeypatch.setattr(datasets_mod.anndata, "read_h5ad", read_h5ad)

    result = loader(sample=sample, save_path="/tmp/save-path")

    assert result is sentinel
    assert read_paths == [f"/tmp/harreman-cache/{fname}"]
    assert calls == [
        {
            "url": f"{BASE_URL}/datasets/{fname}",
            "known_hash": DATASET_HASHES[fname],
            "fname": fname,
            "path": "/tmp/save-path",
            "progressbar": True,
        }
    ]


@pytest.mark.parametrize("species", ["human", "mouse"])
def test_extract_transporter_info_fetches_reachable_files_with_hash(monkeypatch, species):
    calls = _patch_retrieve(monkeypatch, database_mod)
    monkeypatch.setattr(database_mod.pooch, "os_cache", lambda name: f"/cache/{name}")

    def read_csv(path, index_col=0):
        if "Heterodimer" in path:
            return pd.DataFrame(index=["heterodimer"])
        return pd.DataFrame({"Metabolite": ["glucose_imp"], "Gene": ["gene0/gene1"]})

    monkeypatch.setattr(database_mod.pd, "read_csv", read_csv)
    adata = AnnData(
        X=np.ones((2, 2)),
        var=pd.DataFrame(index=["gene0", "gene1"]),
    )

    database_mod.extract_transporter_info(adata, species=species)

    expected_fnames = [
        f"HarremanDB_{species}_extracellular.csv",
        f"HarremanDB_{species}.csv",
        f"Heterodimer_info_{species}.csv",
    ]
    assert calls == [
        {
            "url": f"{BASE_URL}/HarremanDB/{fname}",
            "known_hash": HARREMAN_DB_HASHES[fname],
            "fname": fname,
            "path": "/cache/scvi_harreman",
            "progressbar": False,
        }
        for fname in expected_fnames
    ]


@pytest.mark.parametrize("species", ["human", "mouse"])
def test_extract_lr_pairs_fetches_v2_cellchat_files_with_hash(monkeypatch, species):
    calls = _patch_retrieve(monkeypatch, database_mod)
    monkeypatch.setattr(database_mod.pooch, "os_cache", lambda name: f"/cache/{name}")

    def read_csv(path, index_col=0):
        if "interaction" in path:
            return pd.DataFrame(
                {
                    "ligand": ["gene0"],
                    "receptor": ["gene1"],
                    "annotation": ["secreted signaling"],
                },
                index=["pair0"],
            )
        return pd.DataFrame(index=[])

    monkeypatch.setattr(database_mod.pd, "read_csv", read_csv)
    adata = AnnData(
        X=np.ones((2, 2)),
        var=pd.DataFrame(index=["gene0", "gene1"]),
    )

    database_mod.extract_lr_pairs(adata, species=species)

    expected_fnames = [
        f"interaction_input_CellChatDB_v2_{species}.csv",
        f"complex_input_CellChatDB_v2_{species}.csv",
    ]
    assert calls == [
        {
            "url": f"{BASE_URL}/CellChatDB/{fname}",
            "known_hash": CELLCHAT_HASHES[fname],
            "fname": fname,
            "path": "/cache/scvi_harreman",
            "progressbar": False,
        }
        for fname in expected_fnames
    ]


@pytest.mark.parametrize("species", ["human", "mouse"])
def test_load_metabolic_genes_fetches_from_exampledata_with_hash(monkeypatch, species):
    calls = _patch_retrieve(monkeypatch, autocorr_mod)
    monkeypatch.setattr(autocorr_mod.pooch, "os_cache", lambda name: f"/cache/{name}")
    monkeypatch.setattr(
        autocorr_mod.pd,
        "read_csv",
        lambda path, index_col=0: pd.DataFrame({"0": ["gene0", "gene1"]}),
    )

    genes = autocorr_mod.load_metabolic_genes(species)
    fname = f"{species}_metabolic_genes.csv"

    assert genes == ["gene0", "gene1"]
    assert calls == [
        {
            "url": f"{BASE_URL}/metabolic_genes/{fname}",
            "known_hash": METABOLIC_GENE_HASHES[fname],
            "fname": fname,
            "path": "/cache/scvi_harreman",
            "progressbar": False,
        }
    ]
