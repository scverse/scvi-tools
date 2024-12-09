import os
import subprocess

import pytest
import torch


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
    "label_0",
    batch_key="batch",
)

model = SCANVI(adata, n_latent=10)

model.train(
    max_epochs=100,
    train_size=0.5,
    check_val_every_n_epoch=1,
    accelerator="gpu",
    devices=-1,
    strategy="ddp_find_unused_parameters_true",
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
