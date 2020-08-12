from abc import ABC, abstractmethod


class AbstractModelClass(ABC):
    def __init__(self):
        pass

    @abstractmethod
    def train(self):
        pass

    # @abstractmethod
    # def test(self):
    #     pass
