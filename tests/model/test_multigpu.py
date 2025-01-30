import os
import subprocess

import pytest
import torch

import scvi
from scvi.model import SCANVI, SCVI


@pytest.mark.multigpu
# SCANVI FROM SCVI - reminder: its impossible to debug pytest multigpu work like this
def test_scanvi_from_scvi_multigpu():
    if torch.cuda.is_available():
        adata = scvi.data.synthetic_iid()

        SCVI.setup_anndata(adata)

        datasplitter_kwargs = {}
        datasplitter_kwargs["drop_dataset_tail"] = True
        datasplitter_kwargs["drop_last"] = False

        model = SCVI(adata)

        print("multi GPU SCVI train")
        model.train(
            max_epochs=1,
            check_val_every_n_epoch=1,
            accelerator="gpu",
            devices=-1,
            datasplitter_kwargs=datasplitter_kwargs,
            strategy="ddp_find_unused_parameters_true",
        )
        print("done")
        torch.distributed.destroy_process_group()

        assert model.is_trained
        adata.obsm["scVI"] = model.get_latent_representation()

        datasplitter_kwargs = {}
        datasplitter_kwargs["distributed_sampler"] = True
        datasplitter_kwargs["drop_last"] = False

        print("multi GPU scanvi load from scvi model")
        model_scanvi = scvi.model.SCANVI.from_scvi_model(
            model,
            adata=adata,
            labels_key="labels",
            unlabeled_category="unknown",
        )
        print("done")
        print("multi GPU scanvi train from scvi")
        model_scanvi.train(
            max_epochs=1,
            train_size=0.5,
            check_val_every_n_epoch=1,
            accelerator="gpu",
            devices=-1,
            strategy="ddp_find_unused_parameters_true",
            datasplitter_kwargs=datasplitter_kwargs,
        )
        print("done")
        adata.obsm["scANVI"] = model_scanvi.get_latent_representation()

        torch.distributed.destroy_process_group()

        assert model_scanvi.is_trained


@pytest.mark.multigpu
# SCANVI FROM SCRATCH - reminder: its impossible to debug pytest multigpu work like this
def test_scanvi_from_scratch_multigpu():
    if torch.cuda.is_available():
        adata = scvi.data.synthetic_iid()

        # SCVI.setup_anndata(adata)
        #
        # datasplitter_kwargs = {}
        # datasplitter_kwargs["drop_dataset_tail"] = True
        # datasplitter_kwargs["drop_last"] = False
        #
        # print("multi GPU train")
        # model.train(
        #     max_epochs=1,
        #     check_val_every_n_epoch=1,
        #     accelerator="gpu",
        #     devices=-1,
        #     datasplitter_kwargs=datasplitter_kwargs,
        #     strategy='ddp_find_unused_parameters_true'
        # )
        #
        # torch.distributed.destroy_process_group()
        #
        # assert model.is_trained

        SCANVI.setup_anndata(
            adata,
            "labels",
            "unknown",
            batch_key="batch",
        )

        datasplitter_kwargs = {}
        datasplitter_kwargs["drop_dataset_tail"] = True
        datasplitter_kwargs["drop_last"] = False

        model = SCANVI(adata, n_latent=10)

        print("multi GPU scanvi train from scracth")
        model.train(
            max_epochs=1,
            train_size=0.5,
            check_val_every_n_epoch=1,
            accelerator="gpu",
            devices=-1,
            datasplitter_kwargs=datasplitter_kwargs,
            strategy="ddp_find_unused_parameters_true",
        )
        print("done")

        torch.distributed.destroy_process_group()

        assert model.is_trained


@pytest.mark.multigpu
def test_scvi_train_ddp(save_path: str):
    training_code = """
import torch
import scvi
from scvi.model import SCVI

adata = scvi.data.synthetic_iid()
SCVI.setup_anndata(adata)

model = SCVI(adata)

model.train(
    max_epochs=1,
    check_val_every_n_epoch=1,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true",
)

torch.distributed.destroy_process_group()

assert model.is_trained
"""
    # Define the file path for the temporary script in the current working directory
    temp_file_path = os.path.join(save_path, "train_scvi_ddp_temp.py")

    # Write the training code to the file in the current working directory
    with open(temp_file_path, "w") as temp_file:
        temp_file.write(training_code)
        print(f"Temporary Python file created at: {temp_file_path}")

    def launch_ddp(world_size, temp_file_path):
        # Command to run the script via torchrun
        command = [
            "torchrun",
            "--nproc_per_node=" + str(world_size),  # Specify the number of GPUs
            temp_file_path,  # Your original script
        ]
        # Use subprocess to run the command
        try:
            # Run the command, wait for it to finish & clean up the temporary file
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            os.remove(temp_file_path)
            print(f"Error occurred while running the DDP training: {e}")
            raise
        finally:
            os.remove(temp_file_path)

    launch_ddp(torch.cuda.device_count(), temp_file_path)


@pytest.mark.multigpu
def test_scanvi_train_ddp(save_path: str):
    training_code = """
import torch
import scvi
from scvi.model import SCANVI

adata = scvi.data.synthetic_iid()
SCANVI.setup_anndata(
    adata,
    "labels",
    "unknown",
    batch_key="batch",
)

model = SCANVI(adata, n_latent=10)

datasplitter_kwargs = {}
datasplitter_kwargs["drop_dataset_tail"] = True
datasplitter_kwargs["drop_last"] = False

model.train(
    max_epochs=1,
    train_size=0.5,
    check_val_every_n_epoch=1,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true",
    datasplitter_kwargs=datasplitter_kwargs,
)

torch.distributed.destroy_process_group()

assert model.is_trained
"""
    # Define the file path for the temporary script in the current working directory
    temp_file_path = os.path.join(save_path, "train_scanvi_ddp_temp.py")

    # Write the training code to the file in the current working directory
    with open(temp_file_path, "w") as temp_file:
        temp_file.write(training_code)
        print(f"Temporary Python file created at: {temp_file_path}")

    def launch_ddp(world_size, temp_file_path):
        # Command to run the script via torchrun
        command = [
            "torchrun",
            "--nproc_per_node=" + str(world_size),  # Specify the number of GPUs
            temp_file_path,  # Your original script
        ]
        # Use subprocess to run the command
        try:
            # Run the command, wait for it to finish & clean up the temporary file
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as e:
            os.remove(temp_file_path)
            print(f"Error occurred while running the DDP training: {e}")
            raise
        finally:
            os.remove(temp_file_path)

    launch_ddp(torch.cuda.device_count(), temp_file_path)
