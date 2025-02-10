from __future__ import annotations

from pprint import pprint
from time import time

import bionty as bt
import lamindb as ln
from lamindb.core import datasets

import scvi
from scvi.dataloaders import MappedCollectionDataModule
from scvi.utils import dependencies

#
# # Ingest dataset1
# adata = datasets.small_dataset1(format="anndata")
# curator = ln.Curator.from_anndata(
#     adata,
#     var_index=bt.Gene.symbol,
#     categoricals={
#         "cell_medium": ln.ULabel.name,
#         "cell_type_by_expert": bt.CellType.name,
#         "cell_type_by_model": bt.CellType.name,
#     },
#     organism="human",
# )
# artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
# artifact.features.add_values(adata.uns)
#
# # Ingest dataset2
# adata2 = datasets.small_dataset2(format="anndata")
# curator = ln.Curator.from_anndata(
#     adata2,
#     var_index=bt.Gene.symbol,
#     categoricals={
#         "cell_medium": ln.ULabel.name,
#         "cell_type_by_model": bt.CellType.name,
#     },
#     organism="human",
# )
# artifact2 = curator.save_artifact(key="example_datasets/dataset2.h5ad")
# artifact2.features.add_values(adata2.uns)


@dependencies("lamindb")
def test_lamindb_dataloader_scvi_scanvi(save_path: str = "./test-registries"):
    # import cellxgene_lamin as cxg
    # import pandas as pd
    # this is done once
    # os.system("lamin init --storage ./test-registries --modules bionty")
    # ln.setup.init(name='lamindb_instance_name',storage=save_path)

    # # Create non-curated datasets
    # ln.Artifact(datasets.file_jpg_paradisi05(), key="images/my_image.jpg").save()
    # ln.Artifact(datasets.file_fastq(), key="raw/my_fastq.fastq").save()
    # ln.Artifact.from_df(datasets.df_iris(), key="iris/iris_collection.parquet").save()
    #
    # # Create a more complex case
    # # observation-level metadata
    # ln.Feature(name="cell_medium", dtype="cat[ULabel]").save()
    # ln.Feature(name="sample_note", dtype="str").save()
    # ln.Feature(name="cell_type_by_expert", dtype="cat[bionty.CellType]").save()
    # ln.Feature(name="cell_type_by_model", dtype="cat[bionty.CellType]").save()
    # # dataset-level metadata
    # ln.Feature(name="temperature", dtype="float").save()
    # ln.Feature(name="study", dtype="cat[ULabel]").save()
    # ln.Feature(name="date_of_study", dtype="date").save()
    # ln.Feature(name="study_note", dtype="str").save()
    #
    # ## Permissible values for categoricals
    # ln.ULabel.from_values(["DMSO", "IFNG"], create=True).save()
    # ln.ULabel.from_values(
    #     ["Candidate marker study 1", "Candidate marker study 2"], create=True
    # ).save()
    # bt.CellType.from_values(["B cell", "T cell"], create=True).save()

    # Ingest dataset
    # adata = synthetic_iid()
    # curator = ln.Curator.from_anndata(
    #     adata,
    #     var_index=adata.var_names,
    #     organism="human",
    # )
    # artifact = curator.save_artifact(key="example_datasets/dataset.h5ad")
    # artifact.features.add_values(adata.uns)

    # Ingest dataset1
    adata1 = datasets.small_dataset1(format="anndata")
    curator = ln.Curator.from_anndata(
        adata1,
        var_index=bt.Gene.symbol,
        categoricals={
            "cell_medium": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
        organism="human",
    )
    artifact1 = curator.save_artifact(key="example_datasets/dataset1.h5ad")
    artifact1.features.add_values(adata1.uns)

    # Ingest dataset2
    # adata2 = datasets.small_dataset2(format="anndata")
    # curator = ln.Curator.from_anndata(
    #     adata2,
    #     var_index=bt.Gene.symbol,
    #     categoricals={
    #         "cell_medium": ln.ULabel.name,
    #         "cell_type_by_model": bt.CellType.name,
    #     },
    #     organism="human",
    # )
    # artifact2 = curator.save_artifact(key="example_datasets/dataset2.h5ad")
    # artifact2.features.add_values(adata2.uns)

    # # # a sample dataset
    # df = pd.DataFrame(
    #     {
    #         "CD8A": [1, 2, 3],
    #         "CD4": [3, 4, 5],
    #         "CD14": [5, 6, 7],
    #         "perturbation": ["DMSO", "IFNJ", "DMSO"],
    #     },
    #     index=["sample1", "sample2", "sample3"],
    # )

    # create & save an artifact from a DataFrame -- delete via artifact.delete(permanent=True)
    # artifact = ln.Artifact.from_df(df, description="my RNA-seq").save()

    # describe the artifact
    # artifact.describe()

    # adata1 = synthetic_iid()
    # adata2 = synthetic_iid()
    #
    # adata = ln.core.datasets.anndata_human_immune_cells()
    # #
    # curator = ln.Curator.from_anndata(adata, var_index=adata.var_names)
    # artifact = curator.save_artifact(description="Human immune cells from Conde22")
    # artifact.describe(print_types=True)

    # We must create the Curate object again to ensure that it references the correct AnnData obj
    # curator = cxg.Curator(adata1)
    # artifact = curator.save_artifact(
    #     description=f"dataset curated against cellxgene schema {curator.schema_version}"
    # )
    # artifact.describe()

    # ln.track("NJvdsWWbJlZS0000")
    collection = ln.Collection(artifact1, key="test_collection")
    # adata = collection.load()

    # a test for mapped collection
    # collection = ln.Collection.get(name="covid_normal_lung") #this is from lamin AWS RDS
    # artifacts = collection.artifacts.all()
    # artifacts.df()

    datamodule = MappedCollectionDataModule(collection, batch_size=1024, join="inner")
    model = scvi.model.SCVI(adata=None, registry=datamodule.registry)
    # pprint(model.summary_stats)
    pprint(model.module)
    inference_dataloader = datamodule.inference_dataloader()

    # Using regular adata laoder
    # adata = collection.load()  # try to compare this in regular settings
    # # setup large
    # SCVI.setup_anndata(
    #     adata,
    #     batch_key="assay",
    # )
    # model_reg = SCVI(adata)
    # start_time = time()
    # model_reg.train(max_epochs=10, batch_size=1024)
    # time_reg = time() - start_time
    # print(time_reg)

    start_time = time()
    model.train(max_epochs=1, datamodule=datamodule)
    time_lamin = time() - start_time
    print(time_lamin)

    _ = model.get_elbo(dataloader=inference_dataloader)
    _ = model.get_marginal_ll(dataloader=inference_dataloader)
    _ = model.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model.get_latent_representation(dataloader=inference_dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True)
    model_query = model.load_query_data(
        adata=False, reference_model="lamin_model", registry=datamodule.registry
    )
    model_query.train(max_epochs=1, datamodule=datamodule)
    _ = model_query.get_elbo(dataloader=inference_dataloader)
    _ = model_query.get_marginal_ll(dataloader=inference_dataloader)
    _ = model_query.get_reconstruction_error(dataloader=inference_dataloader)
    _ = model_query.get_latent_representation(dataloader=inference_dataloader)

    adata = collection.load(join="inner")
    model_query_adata = model.load_query_data(adata=adata, reference_model="lamin_model")
    adata = collection.load(join="inner")
    model_query_adata = model.load_query_data(adata)
    model_query_adata.train(max_epochs=1)
    _ = model_query_adata.get_elbo()
    _ = model_query_adata.get_marginal_ll()
    _ = model_query_adata.get_reconstruction_error()
    _ = model_query_adata.get_latent_representation()
    _ = model_query_adata.get_latent_representation(dataloader=inference_dataloader)

    model.save("lamin_model", save_anndata=False, overwrite=True)
    model.load("lamin_model", adata=False)
    model.load_query_data(adata=False, reference_model="lamin_model", registry=datamodule.registry)

    model.load_query_data(adata=adata, reference_model="lamin_model")
    model_adata = model.load("lamin_model", adata=adata)
    model_adata.train(max_epochs=1)
    model_adata.save("lamin_model_anndata", save_anndata=False, overwrite=True)
    model_adata.load("lamin_model_anndata", adata=False)
    model_adata.load_query_data(
        adata=False, reference_model="lamin_model_anndata", registry=datamodule.registry
    )
