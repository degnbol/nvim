#!/usr/bin/env -S kitty +launch
"""Print kitty's text composition gamma and contrast as JSON on stdout.

Imports kitty's internal modules, so it runs only under `kitty +launch`.
"""
import json
import os
import sys

from kitty.config import load_config
from kitty.constants import config_dir, is_macos


def text_composition(strategy: str) -> dict:
    """Gamma and contrast of a `text_composition_strategy` value.

    Mapped and clamped as kitty's `kitty/options/to-c.h` `text_composition_strategy`.

    Args:
        strategy: `platform`, or `gamma [contrast]`.

    Returns:
        gamma: at least 0.01.
        contrast: percent, 0–100.

    Raises:
        ValueError: `legacy`, which kitty composites with another shader path,
            or a malformed value.
    """
    if strategy == "legacy":
        raise ValueError("text_composition_strategy legacy has no gamma and contrast")
    if strategy == "platform":
        gamma, contrast = (1.7, 30.0) if is_macos else (1.0, 0.0)
    else:
        parts = strategy.split(None, 2)
        if not 1 <= len(parts) <= 2:
            raise ValueError(f"text_composition_strategy must be of the form 'number [number]': {strategy}")
        gamma = float(parts[0])
        contrast = float(parts[1]) if len(parts) > 1 else 0.0
    return {"gamma": max(0.01, gamma), "contrast": min(100.0, max(0.0, contrast))}


if __name__ == "__main__":
    opts = load_config(os.path.join(config_dir, "kitty.conf"))
    try:
        print(json.dumps(text_composition(opts.text_composition_strategy)))
    except ValueError as e:
        sys.exit(str(e))
