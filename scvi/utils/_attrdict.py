class attrdict(dict):
    """
    A dictionary that allows for attribute notation (e.g. d.element).

    Based off of https://stackoverflow.com/questions/38034377/object-like-attribute-access-for-nested-dictionary.
    """

    def __init__(self, *args, **kwargs):
        def from_nested_dict(data):
            if not isinstance(data, dict):
                return data
            else:
                return attrdict({key: from_nested_dict(data[key]) for key in data})

        super().__init__(*args, **kwargs)

        for key in self.keys():
            if hasattr(self, key):
                raise ValueError(
                    f"Cannot create attrdict containing key {key} due to conflict with built-in dict attribute."
                )
            self[key] = from_nested_dict(self[key])

        self.__dict__ = self
