from .inference import Inference
from .classifier_inference import ClassifierInference
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
