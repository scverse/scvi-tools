class TunableMeta(type):
    def __getitem__(cls, values):
        if not isinstance(values, tuple):
            values = (values,)
        return type("Tunable_", (Tunable,), dict(__args__=values))


class Tunable(metaclass=TunableMeta):
    pass
