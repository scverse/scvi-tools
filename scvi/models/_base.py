import torch
import os
from abc import ABC, abstractmethod


class AbstractModelClass(ABC):
    def __init__(self, adata, use_cuda):
        assert (
            "scvi_data_registry" in adata.uns.keys()
        ), "Please setup your AnnData with scvi.dataset.setup_anndata(adata) first"
        self.adata = adata
        self.summary_stats = adata.uns["scvi_summary_stats"]
        self.is_trained = False
        self.use_cuda = use_cuda and torch.cuda.is_available()
        self.batch_size = 128
        self.model_summary_string = ""

        self._posterior_class = None
        self._trainer_class = None

    @abstractmethod
    def train(self):
        pass

    def save(self, dir_path):
        # save the model state dict and the trainer state dict only
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        else:
            raise ValueError(
                "{} already exists. Please provide an unexisting directory for saving.".format(
                    dir_path
                )
            )
        torch.save(self.model.state_dict(), os.path.join(dir_path, "model_params.pt"))
        torch.save(
            self.trainer.optimizer.state_dict(),
            os.path.join(dir_path, "optimizer_params.pt"),
        )

    def load(self, dir_path):
        # load state dicts, maybe a class method?
        # Loading scVI model
        model_path = os.path.join(dir_path, "model_params.pt")
        optimizer_path = os.path.join(dir_path, "optimizer_params.pt")
        if self.use_cuda:
            self.model.load_state_dict(torch.load(model_path))
            self.trainer.optimizer.load_state_dict(torch.load(optimizer_path))
            self.model.cuda()
        else:
            device = torch.device("cpu")
            self.model.load_state_dict(torch.load(model_path, map_location=device))
            self.trainer.optimizer.load_state_dict(
                torch.load(optimizer_path, map_location=device)
            )
        self.model.eval()

    def __repr__(self,):
        summary_string = self.model_summary_string + "\nTraining status: {}".format(
            "Trained" if self.is_trained else "Not Trained"
        )
        return summary_string
