#!/usr/bin/env -S kitty +launch
"""Print kitty's cell layout for the configured font as JSON on stdout.

Imports kitty's internal modules, so it runs only under `kitty +launch`.
"""
import json
import os

from kitty import fast_data_types as fdt
from kitty.config import load_config
from kitty.constants import config_dir
from kitty.fonts import FontModification, ModificationUnit
from kitty.fonts.render import set_font_family

# Measure at 10× font_size so kitty's per-size pixel rounding becomes negligible.
SCALE = 10
DPI = 96.0
LAYOUT_METRICS = ("cell_width", "cell_height", "baseline")


def fixed_offsets(modifications: dict[str, FontModification]) -> list[str]:
    """Non-zero pt/px `modify_font` entries for the cell width, cell height or baseline.

    Kitty adds these as a fixed number of pixels, so unlike percent entries they
    change the cell layout relative to the font size.

    Args:
        modifications: kitty's `modify_font` option.

    Returns:
        The entries as written in kitty.conf.
    """
    return [
        repr(m) for m in modifications.values()
        if m.mod_type.name in LAYOUT_METRICS
        and m.mod_value.unit is not ModificationUnit.percent
        and m.mod_value.val != 0
    ]


def cell_fractions() -> dict:
    """Cell layout of the medium face of `<kitty config_dir>/kitty.conf`, after `modify_font`.

    Independent of zoom. Sets kitty's global font options in this process.

    Raises:
        ValueError: `modify_font` shifts the cell width, cell height or baseline by
            a fixed pt/px amount. The layout then depends on zoom.

    Returns:
        Keys below. All values except `aspect` are fractions of the cell height.
        postscript_name: the medium face.
        aspect: cell height / cell width.
        below: distance from the baseline to the cell bottom.
        x: ink height of "x".
        cap: ink height of "H".
        desc: ink height of "p" minus that of "x".
    """
    opts = load_config(os.path.join(config_dir, "kitty.conf"))
    if offsets := fixed_offsets(opts.modify_font):
        raise ValueError(f"zoom-dependent modify_font entries: {', '.join(offsets)}")
    opts = opts._replace(font_size=SCALE * opts.font_size)
    fdt.set_options(opts)
    # Limits as in kitty/fonts/render.py setup_for_testing.
    fdt.sprite_map_set_limits(100000, 100)
    fdt.set_send_sprite_to_gpu(lambda *_: None)
    set_font_family(opts)
    cell_width, cell_height, baseline = fdt.create_test_font_group(opts.font_size, DPI, DPI)
    face = fdt.current_fonts()["medium"]
    ink = {c: face.render_codepoint(ord(c))[2] / cell_height for c in "xHp"}
    return {
        "postscript_name": face.postscript_name(),
        "aspect": cell_height / cell_width,
        "below": (cell_height - baseline) / cell_height,
        "x": ink["x"],
        "cap": ink["H"],
        "desc": ink["p"] - ink["x"],
    }


if __name__ == "__main__":
    print(json.dumps(cell_fractions()))
