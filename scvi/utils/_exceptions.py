from __future__ import annotations

from typing import Any


class InvalidParameterError(Exception):
    """Exception raised for invalid arguments."""

    def __init__(
        self,
        param: str,
        value: Any,
        valid: list | None = None,
        additional_message: str | None = None,
    ):
        self.message = f"Invalid value for `{param}`: {value}."
        if valid is not None:
            self.message += f" Must be one of {valid}."
        if additional_message is not None:
            self.message += f" {additional_message}"
        super().__init__(self.message)

    def __str__(self):
        return self.message
