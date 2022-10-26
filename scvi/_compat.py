try:
    from typing import Literal
except ImportError:
    try:
        from typing import Literal
    except ImportError:

        class LiteralMeta(type):
            """Literal type meta class."""

            def __getitem__(cls, values):
                if not isinstance(values, tuple):
                    values = (values,)
                return type("Literal_", (Literal,), dict(__args__=values))

        class Literal(metaclass=LiteralMeta):
            """Literal type for compatibility."""
