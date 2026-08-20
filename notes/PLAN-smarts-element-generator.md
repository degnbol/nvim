# PLAN: generate the element colours

Every element gets its standard colour, so the palette covers the whole table
rather than nine elements: too many to type, and the standard values are not
usable as published. Both facts are measured below, and both are what the
generator exists to handle.

How the colours reach the screen is `notes/PLAN-smarts-element-perf.md`, which
supersedes this note's query-per-element design.

## Why it has to be computed

**Size.** The grammar's bracket alphabet holds 117 symbols, of which 91 are
metals — counting polonium, astatine and tennessine as non-metals, which is a
convention rather than a fact, so state it in the generator and stop arguing.
Jmol's table stops at meitnerium, leaving `Ds`, `Rg`, `Cn`, `Fl`, `Mc`, `Lv`,
`Ts` and `Og` with no colour; all eight are in the suspect set below, so the
lookup stays total either way. 109 symbols have a Jmol colour.

**Readability.** Raw Jmol values are unusable on nearly half the table. Against
the two generated backgrounds, taking 2.3:1 as a trial floor:

| set | below the floor on light | on dark | on either |
| --- | --- | --- | --- |
| all coloured | 48 | 8 | 56 |
| metals only | 37 | 5 | 42 |

Praseodymium `#d9ffc7` reads 1.01:1 on `#faf4ed`, i.e. invisible; francium
`#420066` reads 1.07:1 on `#152528`. The existing palette already handles this by
hand for two elements ("sulfur and fluorine are darkened…"). At 42 elements the
adjustment has to be a rule, and a rule is a script.

**The floor is 2.0, per background.** Measured on the shipped palette, the
weakest colour is fluorine at 2.24:1 on light and iodine at 2.02:1 on dark —
phosphorus (2.31 light) is not the weakest, as an earlier draft of this plan
claimed. If the floor means "no worse than what already ships", it is 2.0, and it
differs by background.

## Source of truth

- **Colours**: `ase.data.colors.jmol_colors`, indexed by atomic number, with
  `chemical_symbols` and `atomic_names` from the same package — verified entry for
  entry against Jmol's published table, no mismatches. A PEP 723 uv script, so
  there is no environment to maintain.
- **Which symbols exist**: `modules/tree-sitter-smarts/symbols.js`, parsed rather
  than re-listed. It is itself generated from RDKit probes, so the palette cannot
  drift from what the grammar can produce — including the four disjoint alphabets,
  which decide whether an element has a legal aromatic spelling — the one thing
  the text lookup needs beyond the symbol itself.

## The lightness rule

Per background, in CIELAB, keeping hue and chroma:

1. If the Jmol colour already clears the floor, **keep it unchanged** — so iron
   stays `#e06633` and nitrogen `#3050f8`, and the palette is recognisably CPK.
2. Otherwise move lightness **away from the background** — down for the light
   theme, up for the dark one. Direction is set by the background, not by the
   colour's own luminance: getting that backwards pushes a near-white toward
   white and makes it worse (measured, on the first attempt at this rule).
3. Compress the whole unreadable tail into a 26-point band just inside the
   boundary, preserving order, rather than clamping each element to the boundary
   independently. Independent clamping is what collapses a block of colours onto
   one value.
4. Reduce chroma only where lightness alone cannot reach the floor.

Measured over all 109 elements with steps 1–3: on light, minimum contrast 2.03
and seven pairs below the 2.3 JND; on dark, two pairs. The sub-JND pairs are
`Ca`/`Sr` and the `Eu`/`Gd`/`Tb`/`Dy` cluster — pairs Jmol itself renders nearly
identically, so the rule is not what makes them indistinguishable. Step 4 is
load-bearing on dark, where deep purples (`Fr` `#420066`, `Cs` `#57178f`) cannot
reach 2.0 at full chroma: without it the dark column bottoms out at 1.23.

**The achromatic block survives this rule**, which is the thing to check, since
Jmol distinguishes hydrogen, scandium, titanium, vanadium, silver, platinum and
carbon by lightness alone — zero chroma, no hue to fall back on. Under step 3
they land on `#6c6c6c`, `#808080`, `#9c9fa4`, `#a6a6ab`, `#a0a0a0`, `#90909f` and
`#909090`: compressed, ordered, and mostly still a JND apart. Under independent
clamping they collapse onto a single grey.

**Two values per element, one per background.** Not because a single value cannot
clear the floor on both — at 2.0 it can, and even 3:1 leaves a usable luminance
window — but because it squeezes all 109 elements into L\* 38–67, where the
per-background columns get L\* 16–71 and 21–100. `vim.o.background` is set by the
colorscheme before `ColorScheme` fires (verified), and `hi.onColorScheme` already
re-runs there, so the palette picks a column at no cost.

## Elements that cannot occur in a structure

Some symbols the parsers accept name elements no one has ever made a compound of:
astatine and francium, whose longest-lived isotopes last hours and minutes, and
everything from fermium up, which exists in tracer amounts and single atoms.
Twenty symbols in the grammar's alphabet: `At`, `Fr`, `Fm`, `Md`, `No`, `Lr`,
`Rf`, `Db`, `Sg`, `Bh`, `Hs`, `Mt`, `Ds`, `Rg`, `Cn`, `Fl`, `Mc`, `Lv`, `Ts`,
`Og`.

A structure holding one is a typo or a joke, and that is worth saying out loud
rather than colouring it like any other element — `@constant` is what every
identity already takes, so it reads as "nothing to see here". One
`@chem.suspect`, linked to `DiagnosticWarn`, says the opposite.

This also disposes of the eight elements Jmol leaves uncoloured: they fall inside
this set, so the colour table needs no fallback. Note the set is larger than
Jmol's gap and cuts across it — `Fr` and `At` both have CPK colours they will not
be using.

Two of the twenty are already spelled out in
`after/queries/smarts/highlights.scm`, where astatine and tennessine sit outside
the halogen groups on the grounds that they have no place in a real structure.
That comment becomes this group's justification instead of a local exception.

No library encodes "has an isolable compound", so this is the one hand-maintained
list in the generator. It is small, it is stable, and the alternative boundaries
(past uranium, past lawrencium) are defensible too — the list is the place to
argue about it.

## Output

One generated file, carrying a "generated — do not edit" header like
`symbols.js`, and committed so a machine without uv still gets colours. One row
per element, plus the suspect symbols and which symbols have a legal aromatic
spelling:

```lua
{ symbol = "Fe", z = 26, name = "iron", light = "#…", dark = "#…" },
```

The group is `@chem.element.<name>` rather than `@chem.<name>`, so the documented
`@`-group fallback (`:h treesitter-highlight-groups`) degrades an element to one
defined colour if the palette has not loaded, and element identity is named apart
from the style captures `@chem.bond`, `@chem.anchor`, `@chem.reaction`.

**Hydrogen is the generator's one exception.** It has no symbol node: RDKit reads
`[H]` as the hydrogen atom (round-tripping it to `[#1]`, and `[H+]` to `[#1&+]`)
and as a count of attached hydrogens only where a digit follows it (`[H1]`) or
another primitive precedes it (`[nH]` → `[n&H1]`). No text predicate can tell
those apart, so the atom is matched structurally in
`after/queries/smarts/highlights.scm` — three shapes, keyed on the anchors that
exclude a `(count)` sibling and a preceding element. Hydrogen keeps its colour
row, and only its `#1` spelling is looked up by text.

## The palette's trigger

The groups are defined from an ftplugin today, and a fence has no ftplugin.
Measured: in a markdown buffer holding a ```smiles fence, the element captures do
reach the fence but `@chem.carbon` and `@chem.bond` are undefined, so a fenced
structure highlights with no element colours at all.

The palette belongs with the parser wiring, in `plugin/treesitter.lua` beside
`language.register("smarts", "smiles")` — **not** in the `User TSUpdate` block of
`lua/plugins/treesitter.lua`, which an earlier draft of this plan proposed. That
event is emitted by nvim-treesitter's `install()`, and the call is behind
`if vim.fn.executable("tree-sitter") == 1`, so on a machine without the CLI it
fires zero times and the palette would never be defined at all — worse than
today. Fixing the trigger is a prerequisite: generating a hundred groups behind a
hook half the consumers never reach would multiply the bug.

## Costs — wrong conclusion, kept as the record

Per-line `iter_captures` over a 50-line window, the shape the highlighter uses:

| query | organic input | bracket-heavy input |
| --- | --- | --- |
| today's 9 patterns | 2.1 ms | — |
| ~109 patterns, per-element node types | 10.3 ms | 45.9 ms |
| ~109 patterns, one uniform 4-type union | 69.3 ms | 49.9 ms |

This measured the wrong axis and concluded the per-element node-type lists were
load-bearing. They only ever helped `organic_symbol`; every bracket atom is
`element_symbol` or `atomic_number` and so appears in all ~97 patterns either
way, which is why the bracket column barely moved. A pattern per element is
O(elements) per atom by construction — see
`notes/PLAN-smarts-element-perf.md`, which drops the query layer entirely.

The ~109 `nvim_set_hl` calls cost 0.09 ms, and `:hi` grows by that many entries
in a private namespace no colorscheme reads.

## Tests

- `tests/plenary/smarts_spec.lua` — `[Fe]` and `[#26]` take the same group; `Sn`
  outside brackets is sulfur + nitrogen and takes neither iron nor tin (`Fe` is
  *not* the example for this: bare `Fe` lexes as `F` plus an ERROR, and on line 1
  as a header field); the groups are defined in a markdown fence as well as in a
  `.smi` buffer. Invariant to assert: exactly one `@chem.element.*` group per
  glyph — every element glyph also carries the grammar's own `@constant`, so
  "one group per glyph" is the wrong assertion.
- `tests/plenary/colors_spec.lua` — the light/dark switch, per `tests/README.md`,
  which puts highlight cases there rather than in a filetype spec.
- A data invariant over the generated table: every colour clears its own
  background's floor, and no two elements are within ΔE 2. A cheap loop, and it
  is what catches a regression in the lightness rule.
- A drift check that the generated element set still matches `symbols.js`.

## Open questions

1. **Chlorine** ships CDK `#1f7f1f` at 4.67:1 on light; the Jmol-derived value is
   `#29bc22` at 2.31:1. Accept the downgrade for one source of truth, or keep
   hand-overrides for the organic nine?
2. Where the `@chem.suspect` boundary sits: no isolable compound (the 20 symbols
   above), past lawrencium, or past uranium. The first is the chemical fact; the
   others are rounder numbers.
