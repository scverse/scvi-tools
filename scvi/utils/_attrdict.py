from ml_collections.config_dict import FrozenConfigDict


class attrdict(FrozenConfigDict):
    """
    Dummy class to allow for backwards compatibility with old attrdict
    Now uses ml-collections FrozenConfigDict.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
