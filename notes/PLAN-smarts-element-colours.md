# PLAN: a colour for metals in chemical line notation

The element layer now colours carbon, nitrogen, oxygen, sulfur, phosphorus and
the four halogens. Everything else takes `@constant`, the capture element
identities share with bond orders, so `[Fe]`, `[Pt]` and `[Na]` read as the same
neutral colour a `=` does — and a metal centre is the one thing in a structure a
reader most wants to spot.

## There is no standard metal colour (measured)

RDKit ships three atom palettes, readable with `MolDrawOptions.getAtomPalette()`:

| palette | element entries | covers |
| --- | --- | --- |
| default | 10 | H C N O F P S Cl Br I |
| CDK | 11 | the same plus B |
| Avalon | 10 | the same set, different hues |

Every palette holds three non-element keys: `-1` the fallback for any element it
does not list, `0` the dummy-atom `*` colour, and `201`. **The fallback is
`#000000`, and so is carbon** — RDKit does not draw metals in a shared "other"
colour, it draws them as carbon. So the reference implementation argues for
`@chem.metal` linking to `Normal`, which is what already happens by default.
Any metal colour here is a house choice, not a convention.

Jmol's per-element CPK table (<https://jmol.sourceforge.net/jscolors/>) does
assign every metal a hue, and those hues land on what is already allocated: iron
`#e06633` against phosphorus `#ff8000`, copper `#c88033` beside it, magnesium
`#8aff00` and calcium `#3dff00` against the halogen greens.

## Two options

**A. One `@chem.metal`.** Says "a metal centre is here" and nothing more, which
is the useful signal at this palette's resolution, and needs one new colour.

**B. Jmol CPK for the metals that turn up in organometallic chemistry** — Li Na
K Mg Al Fe Cu Zn Pd Pt Sn, plus the metalloids B Si — with the remaining ~80
staying `@constant`. Standard colours throughout, at the cost of overlapping
hues (iron beside phosphorus) and a per-element list to maintain.

## Query shape

Exclusion is only sound for `(element_symbol)`, whose alphabet is a closed
117-symbol set (`symbols.js`):

```scheme
((element_symbol) @chem.metal
 (#not-any-of? @chem.metal
  "H" "He" "B" "C" "N" "O" "F" "Ne" "Si" "P" "S" "Cl" "Ar" "Ge" "As" "Se" "Br"
  "Kr" "Sb" "Te" "I" "Xe" "At" "Rn" "Ts" "Og")
 (#set! priority 101))
```

That list is exactly {non-metals ∪ metalloids ∪ noble gases}, which no file in
the module states — say so in the query comment, and note that polonium falls
through to metal while the other metalloids do not.

**`(atomic_number)` needs an enumeration instead.** Its text is unbounded, not a
closed set: `[#0]`, `[#06]` and `[#200]` all parse as `atomic_number` (checked
with `tree-sitter parse`), and RDKit accepts all three — `[#0]` is the dummy-atom
query. Under `#not-any-of?` every one of them becomes a metal. So the digit
spellings want `#any-of?` with the 92 metal numbers.

Two details, neither of them a collision:

- `"H"` in the list above is unreachable — `[H]` parses as `hydrogen`, never
  `element_symbol` — but `"#1"` is still needed.
- `[!Fe]` colours `Fe` as a metal, since negation wraps the primitive and leaves
  the symbol node intact. Same as the existing `[!C]` behaviour.
- `Nh` is absent from the grammar's alphabet (SMARTS reads `[Nh]` as N + h-count),
  so `#113` has no symbol spelling.

## Contrast

No single hex clears 4:1 on both generated backgrounds (`#faf4ed`, `#152528`),
and the existing element colours do not either: nitrogen is 2.73 on dark,
phosphorus 2.31 on light. So a metal colour either accepts the same weakness on
one side, or `hi.onColorScheme` branches on `vim.o.background` and gives two
values — the hook already re-runs per colorscheme.

Pick the value by measuring, not by name: a steel-blue `#4682b4` is ΔE 24.6 from
the light theme's `Normal` (`#575279`) and 26.5 from its `Comment` (`#5f5695`),
both blue-violets, where every existing element-vs-carbon pair is ≥68.

## Where it goes

- `after/queries/smarts/highlights.scm` — the new patterns, beside the nine
  element ones. Element identity is node text rather than node type, so this
  stays a query-layer concern.
- `lua/chem_palette.lua` — the group, and the header comment's note on which
  values are darkened and why.

## Reproducibility

The palette numbers are read out of the installed RDKit rather than transcribed:

```python
from rdkit.Chem.Draw import rdMolDraw2D
opts = rdMolDraw2D.MolDrawOptions()
opts.useCDKAtomPalette()
{z: "#%02x%02x%02x" % tuple(round(c * 255) for c in tuple(rgb)[:3])
 for z, rgb in opts.getAtomPalette().items()}
```

## Open questions

1. Option A or B.
2. If A, the hex — and whether it differs by background.
3. Metalloids: with the metals, or left as `@constant`.
</content>
</invoke>
