"""Shared utility functions for Harreman."""

from __future__ import annotations

import torch


def _resolve_device(device: torch.device | str) -> torch.device:
    """Resolve device string and fall back to CPU if CUDA is unavailable."""
    if isinstance(device, str) and device == "auto":
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if isinstance(device, str):
        device = torch.device(device)
    if device.type == "cuda":
        try:
            torch.tensor([0.0], device=device)
        except RuntimeError:
            device = torch.device("cpu")
    return device
