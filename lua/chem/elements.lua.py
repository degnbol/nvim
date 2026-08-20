#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.11"
# dependencies = ["ase", "coloraide"]
# ///
"""A colour for every element a chemical structure spells.

Reads `symbols.js` for which symbols the grammar accepts and which of them have a
legal aromatic (lowercase) spelling, and writes `elements.lua`.

Jmol's published CPK colours are unreadable on nearly half the table —
praseodymium `#d9ffc7` reads 1.01:1 on the light background, francium `#420066`
1.07:1 on the dark one — so each column keeps the published colour where it
already clears the floor and corrects its lightness where it does not. The two
columns are separate because one value legible on both backgrounds would fit
every element into a third of the lightness range.
"""

import json
import re
from pathlib import Path

from ase.data import atomic_names, atomic_numbers
from ase.data.colors import jmol_colors
from coloraide import Color

CONFIG = Path(__file__).parents[2]

# `Normal`'s background in the light and dark colorschemes this config ships
# (dawnfox and terafox, hardlinked in as `colors/generated.lua`). Change a
# colorscheme and the whole palette has to be regenerated against the new pair.
BACKGROUNDS = {"light": Color("#faf4ed"), "dark": Color("#152528")}

# Weakest contrast the hand-written palette already shipped: fluorine at 2.24:1
# on light, iodine at 2.02:1 on dark. Anything above this would move colours
# that have been readable all along.
FLOOR = 2.0

# L* the corrected colours are spread over, starting at the boundary and running
# away from the background. Wide enough to keep the seven achromatic elements —
# the ones Jmol tells apart by lightness alone — about a JND apart, and narrow
# enough that the far end stays off black and white.
BAND = 26.0

# Elements no one has ever made an isolable compound of: astatine and francium,
# whose longest-lived isotopes last hours and minutes, and everything from
# fermium up, which exists in tracer amounts and single atoms. A structure
# holding one is a typo or a joke, so they take a warning rather than a colour.
# No library encodes this, and the neighbouring boundaries (past uranium, past
# lawrencium) are defensible too — argue about it here.
SUSPECT = ("At", "Fr", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs",
           "Mt", "Ds", "Rg", "Cn", "Fl", "Mc", "Lv", "Ts", "Og")


def alphabets(source: str) -> dict[str, list[str]]:
    """The symbol sets a generated JavaScript module declares.

    Args:
        source: module text, one `const name = [...];` per set.

    Returns:
        Set name → its members.
    """
    return {name: json.loads(members)
            for name, members in re.findall(r"const (\w+) = (\[[^\]]*\]);", source)}


def boundary_lightness(background: Color) -> float:
    """The lightness at which a colour first clears the floor against a background.

    WCAG relative luminance is CIE Y and CIE L* is a function of Y alone, so one
    L* bounds the readable range whatever the hue and chroma — which is also why
    correcting a colour means moving its lightness and nothing else.

    Args:
        background: the colour read against.

    Returns:
        Boundary L*, in CIE Lab D65.
    """
    contrasted = (background.luminance() + 0.05)
    luminance = (contrasted / FLOOR if is_light(background)
                 else contrasted * FLOOR) - 0.05
    return Color("xyz-d65", [0, luminance, 0]).convert("lab-d65")["lightness"]


def is_light(background: Color) -> bool:
    """Whether a colour is light, and so needs darker text on it.

    Args:
        background: the colour to judge.

    Returns:
        True if it sits in the lighter half of the lightness range.
    """
    return background.convert("lab-d65")["lightness"] > 50


def rounded(colour: Color, background: Color) -> str:
    """Round to 8 bits per channel without dropping below the floor.

    Args:
        colour: the colour to round, in any space.
        background: the colour it is read against.

    Returns:
        A `#rrggbb` string clearing the floor against the background.
    """
    away = -1 if is_light(background) else 1
    lch = colour.convert("lch-d65")
    while True:
        hexed = lch.convert("srgb").fit(method="lch-chroma").to_string(hex=True)
        if Color(hexed).contrast(background, method="wcag21") >= FLOOR:
            return hexed
        lch["lightness"] += away * 0.2


def readable(colours: dict[str, Color], background: Color) -> dict[str, str]:
    """Correct a palette's illegible colours against one background.

    A colour already clearing the floor is returned unchanged, so iron stays
    `#e06633` and the palette is recognisably CPK. The rest keep their hue and
    are reflected through the boundary lightness into a band just inside it: the
    closer a colour started to the boundary the less it moves, so the corrected
    tail joins the untouched colours continuously instead of collapsing onto the
    boundary. Reflecting spreads them over the band rather than clamping, which
    is what keeps the achromatic elements — distinguished by lightness alone —
    apart. Chroma the new lightness puts outside sRGB is dropped by gamut
    mapping, without which the deep purples cannot move at all.

    Args:
        colours: key → colour.
        background: the colour they are read against.

    Returns:
        key → `#rrggbb` string.
    """
    away = -1 if is_light(background) else 1
    boundary = boundary_lightness(background)
    excess = {key: max(away * (boundary - colour.convert("lch-d65")["lightness"]), 0.0)
              for key, colour in colours.items()}
    worst = max(excess.values())

    corrected = {}
    for key, colour in colours.items():
        if not excess[key]:
            corrected[key] = colour.to_string(hex=True)
            continue
        lch = colour.convert("lch-d65")
        lch["lightness"] = boundary + away * BAND * excess[key] / worst
        corrected[key] = rounded(lch, background)
    return corrected


def row(symbol: str, columns: dict[str, dict[str, str]], aromatic: bool) -> str:
    """One `ChemElement` table constructor.

    Args:
        symbol: an element symbol from the bracket alphabet.
        columns: background name → symbol → `#rrggbb`.
        aromatic: whether the symbol's lowercase spelling parses.

    Returns:
        An indented Lua table constructor, trailing comma included.
    """
    fields = [f'symbol = "{symbol}"',
              f"z = {atomic_numbers[symbol]}",
              f'name = "{atomic_names[atomic_numbers[symbol]].lower()}"',
              *(f'{name} = "{column[symbol]}"' for name, column in columns.items())]
    if aromatic:
        fields.append("aromatic = true")
    return "    { " + ", ".join(fields) + " },"


sets = alphabets((CONFIG / "modules/tree-sitter-smarts/symbols.js").read_text())
symbols = [s for s in sets["bracketElement"] if s not in SUSPECT]
aromatics = set(sets["bareAromatic"] + sets["bracketAromatic"])
jmol = {symbol: Color("srgb", list(jmol_colors[atomic_numbers[symbol]])) for symbol in symbols}
columns = {name: readable(jmol, background) for name, background in BACKGROUNDS.items()}

rows = [row(symbol, columns, symbol.lower() in aromatics) for symbol in symbols]
suspects = [f"        {symbol} = {atomic_numbers[symbol]}," for symbol in SUSPECT]
(CONFIG / "lua/chem/elements.lua").write_text("\n".join([
    f"-- Generated by {Path(__file__).name} — do not edit.",
    "-- Jmol CPK colours, corrected for contrast against each background.",
    "--- @class ChemElement",
    "--- @field symbol string",
    "--- @field z integer",
    "--- @field name string",
    "--- @field light string",
    "--- @field dark string",
    "--- @field aromatic? boolean whether the lowercase spelling parses",
    "",
    "--- @class ChemElements",
    "--- @field [integer] ChemElement",
    "--- @field suspect table<string, integer> symbol -> atomic number",
    "",
    "--- @type ChemElements",
    "return {", *rows, "",
    "    -- Elements with no isolable compound, listed apart because they have no",
    "    -- colour: one identity among others would say there is nothing to see.",
    "    suspect = {", *suspects, "    },",
    "}", ""]))
