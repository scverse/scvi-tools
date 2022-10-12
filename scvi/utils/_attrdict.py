from ml_collections.config_dict import FrozenConfigDict


class attrdict(FrozenConfigDict):
    """
    A dictionary that allows for attribute-style access.

    Dummy class that allowsfor backwards compatibility with the previous custom ``attrdict``.
    Inherits from :class:`~ml_collections.config_dict.FrozenConfigDict`.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
