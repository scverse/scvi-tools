from .inference import Inference
from .classifier_inference import ClassifierInference
from .variational_inference import (
    VariationalInference,
    AlternateSemiSupervisedVariationalInference,
    JointSemiSupervisedVariationalInference,
    VariationalInferenceFish
)
from .experimental_inference import adversarial_wrapper, mmd_wrapper

__all__ = ['Inference',
           'ClassifierInference',
           'VariationalInference',
           'AlternateSemiSupervisedVariationalInference',
           'JointSemiSupervisedVariationalInference',
           'VariationalInferenceFish',
           'adversarial_wrapper',
           'mmd_wrapper']
