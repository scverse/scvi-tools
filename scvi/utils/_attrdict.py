from copy import deepcopy


class attrdict(dict):
    """
    A dictionary that allows for attribute notation (e.g. d.element).

    Parameters
    ----------
    recursive
        If True, recursively converts nested dictionaries into :class:`~scvi.utils.attrdict` objects.

    Notes
    -----
    Based off of https://stackoverflow.com/questions/38034377/object-like-attribute-access-for-nested-dictionary.
    """

    def __init__(self, *args, recursive: bool = False, **kwargs):
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
            if recursive:
                self[key] = from_nested_dict(self[key])
            else:
                self[key] = deepcopy(self[key])

        self.__dict__ = self

    def __repr__(self) -> str:
        return f"attrdict({super().__repr__()})"
