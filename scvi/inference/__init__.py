from .classifier_inference import ClassifierInference
from .inference import Inference
from .variational_inference import (
    VariationalInference,
    AlternateSemiSupervisedVariationalInference,
    JointSemiSupervisedVariationalInference
)

__all__ = ['Inference',
           'ClassifierInference',
           'VariationalInference',
           'AlternateSemiSupervisedVariationalInference',
           'JointSemiSupervisedVariationalInference']
