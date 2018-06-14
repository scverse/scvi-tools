''' Defines the abstract types and required signatures for models '''

from abc import ABC, abstractmethod


class BaseModel(ABC):
    @abstractmethod
    def forward(self, x, local_l_mean, local_l_var, batch_index=None, y=None):
        '''
        :return: reconst_loss, kl_divergence
        '''
        pass

    @abstractmethod
    def get_latents(self, x, y=None):
        '''
        :param x: sample_batch
        :param y: labels
        :return: an array [z1, z2] the sampled variables from bottom to top
        '''
        pass


class SemiSupervisedModel(BaseModel):
    def classify(self, x):
        '''
        :return: a softmax probability over all labels for every sample in x.
        '''
        pass
