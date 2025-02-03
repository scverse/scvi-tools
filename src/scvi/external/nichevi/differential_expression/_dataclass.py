from dataclasses import dataclass, field
from typing import Union

import numpy as np
import pandas as pd
from rich import print
from sklearn.gaussian_process import GaussianProcessClassifier


@dataclass
class DifferentialExpressionResults:
    """Dataclass for storing the results of the differential expression analysis,
    including the GP classifier
    """  # noqa: D205

    gpc: GaussianProcessClassifier
    g1_g2: pd.DataFrame
    g1_n1: pd.DataFrame
    n1_g2: pd.DataFrame
    n1_n2: Union[pd.DataFrame, None] = field(default=None)  # noqa: UP007
    n1_index: Union[np.array, None] = field(default=None)  # noqa: UP007
    n2_index: Union[np.array, None] = field(default=None)  # noqa: UP007

    def gpc_info(self):
        """Print the log marginal likelihood value and the kernel
        of the Gaussian Process Classifier

        """  # noqa: D205
        print("Training score: ", self.gpc.train_score_)
        print("Marginal likelihood: ", self.gpc.log_marginal_likelihood_value_)
        print("Kernel: ", self.gpc.kernel_)
