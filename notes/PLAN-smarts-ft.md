# PLAN: chemical line notation filetypes and grammar

Highlight SMILES, SMARTS, SMIRKS, CXSMILES and their reaction forms — in
dedicated files, and where they are embedded in other languages.

## Scope

One tree-sitter grammar named **`smarts`**. Every member of the family is either
a subset of SMARTS or SMARTS plus a layer the grammar already carries: SMILES is
very nearly a syntactic subset (four exceptions, measured below), SMIRKS is a
subset of reaction SMARTS, and CXSMARTS is SMARTS with the same CX layer as
CXSMILES. So `smarts` names the grammar's actual reach; `smiles` would name a
subset.

Two filetypes, `smiles` and `smarts`, share that grammar. A `.smarts` file is not
written in SMILES, and the filetype is what records which dialect a file is in.
Each additional filetype is one registration line and no grammar work.

RDKit is the arbiter of validity and of tokenisation: it is the reference
implementation these files are read by, so where it diverges from Daylight or
OpenSMILES the grammar follows RDKit. Concretely, `(` always follows an atom —
RDKit rejects the top-level component grouping `(C).(C)` that the Daylight SMARTS
spec allows.

RDKit reaches the family through two mutually exclusive entry points, and the
grammar spans both: `MolFromSmarts` takes a single structure and rejects `>`,
`ReactionFromSmarts` requires `>` and rejects a lone structure. Both accept the
CX layer.

The grammar accepts the union, so a SMARTS-only construct in a `.smi` file
highlights cleanly. Ring-closure pairing, valence and aromaticity perception are
semantic and stay with RDKit.

## Language behaviour (measured)

From `rdkit` probes, parse-only (`sanitize=False`) where chemistry validity would
otherwise mask a lexing result.

**SMARTS is a syntactic superset of SMILES with four exceptions.** Probing all
118 element symbols in brackets through both parsers:

```
[Mc] [Ts] [Og]   accepted as SMILES, rejected by MolFromSmarts
[Nh]             SMILES: nihonium.  SMARTS: [N&h] — nitrogen + implicit-H count
```

Everything else agrees, including the constructs that look like collisions: `C$C`
(quadruple bond, while `$(` opens recursive SMARTS — no clash, `$(` occurs only
inside brackets), the extended chirality classes `[C@TH1]` `[C@AL1]` `[S@TB5]`
`[Co@OH12]`, `%(100)` ring closures, and atom maps.

**Take the SMARTS reading of `Nh`.** `h` is everyday SMARTS; nihonium is not.
`[Nh]` and `[nH]` differ in case, so this costs only nihonium. Record it in the
corpus as a documented divergence.

**SMIRKS needs no separate handling.** Per the RDKit maintainer, "the RDKit
Reaction SMARTS parser should parse, and probably correctly interpret, any valid
SMIRKS" ([discussion 5168](https://github.com/rdkit/rdkit/discussions/5168)). The
divergence there is semantic — SMIRKS carries explicit hydrogens, which affects
matching against H-suppressed molecules, not tokenisation.

**The CX layer applies to SMARTS, not just SMILES.** It is parsed, not skipped:

```
MolFromSmarts('[#6][#6][#8] |$foo;;$|')  -> atomLabel ['foo', '-', '-']
MolFromSmarts('CCO |(0,0,;1,0,;2,0,)|')  -> conformers=1
ReactionFromSmarts('[C:1][O:2]>>[C:1][O:2] |$a;b;c;d$|')  -> ok
```

So CXSMARTS is covered by the one `cx_layer` rule, on molecules and reactions
alike.

**Inside brackets, element symbols win by longest match** — in the SMARTS parser
too, not just SMILES:

```
[Ra] -> [#88]  radium, not R&a      [Sc] -> [#21]
[As] -> [#33]  arsenic, not A&s     [Ar] -> [#18]
```

Tree-sitter's longest-match lexer reproduces this for free (`Nh` excepted).

**Outside brackets the alphabet is a different, much smaller set.** This is the
fact the grammar rests on:

```
bare aliphatic:  B C N O F P S Cl Br I
bare aromatic:   b c n o p s
bracket-only:    the other 108, plus aromatic si as se te
```

118 symbols, 104 of them two letters, none longer. A bare two-letter sequence is
therefore **two atoms**, not an element:

```
Sn         ->  ['S', 'N']   sulfur + nitrogen, not tin
[Sn]       ->  Sn           tin
Sc1ccccc1  ->  C6H6S        thiophenol: S + aromatic c, not scandium
Se, Si, Ba ->  NULL         bracket-only
```

The same two characters lex as one token or two depending on bracket depth, so
the element tables must stay disjoint — sharing them produces wrong highlights,
not merely permissive ones.

**Record layout.** Space and tab are interchangeable as the separator, the CX
layer sits between structure and name, and the name is the whole remainder
including internal spaces:

```
'CCO ethanol'            -> name 'ethanol'
'CCO\tethyl alcohol'     -> name 'ethyl alcohol'
'CCO |$_R1;;$| ethanol'  -> CX '|$_R1;;$|', name 'ethanol'
```

The field terminator is **any whitespace**, not a tab.

**A header line is the default.** `SmilesMolSupplier`'s signature is
`titleLine=True`, so the first line of a `.smi` file is expected to be column
names, not a structure. The grammar needs a rule for it.

**Filetype override.** `.smi` resolves to `mib` built-in. An `extension` entry in
`vim.filetype.add` overrides it with no `priority` needed (verified under
`--clean`); the priority caveat applies to `pattern` functions, not extensions.
This drops the vim-regex `mib` highlighting that `.smi` currently gets, including
on a machine with no tree-sitter CLI.

## Comments

`#` opens a comment **only at the start of a line**, leading whitespace allowed.
`C#C` is ethyne, so `#` mid-line is a triple bond; at line start it is
unambiguous because no structure may begin with a bond. Text after the structure
is already the name field.

## Grammar shape

**`extras: []`.** The record layer is whitespace-significant, so no implicit
whitespace skipping — every separator and newline is threaded explicitly. This is
the structural decision everything below depends on.

Outer layer, the `.smi` file convention:

```
source_file → header? line*
line        → comment | record | blank
blank       → ws? eol
comment     → ws? '#' /[^\r\n]*/ eol
record      → structure (ws cx_layer)? (ws name)? ws? eol
header      → header_field (sep header_field)+ eol
cx_layer    → token('|' /[^|\r\n]*/ '|')
name        → /[^|\s][^\r\n]*/
eol         → /\r?\n/
```

- `cx_layer` must be a single token, otherwise `name` outlexes the `'|'` by
  longest match and the CX layer is silently absorbed with no ERROR node.
- `name` may not begin with `|` for the same reason.
- A first line is ambiguous between `header` and `record` when it parses as both
  (`CCO`); give `record` the higher dynamic precedence so structures win and only
  a non-structure first line is a header. Confirm with `tree-sitter generate`.
- `\r?\n` throughout rather than requiring `fileformat=unix`.

Inner layer, OpenSMILES with the SMARTS extensions folded in:

```
structure     → side | side '>' side | side '>' side '>' side
side          → chain | ε at the use site, never as a rule
chain         → branched_atom (bond? chain | '.' chain)?
branched_atom → atom ring_bond* branch*
branch        → '(' bond? chain ')'
atom          → organic_symbol | aromatic_symbol | '*' | 'A' | 'a' | bracket_atom
ring_bond     → bond? (digit | '%' digit digit | '%(' digit+ ')')
bond          → bond_prim (('&' | ',' | ';')? bond_prim)*
bond_prim     → '!'? one-of(- = # $ : / \ ~ @)
bracket_atom  → '[' atom_expr ']'
atom_expr     → and_low
and_low       → or_expr (';' or_expr)*        -- lowest precedence
or_expr       → and_high (',' and_high)*
and_high      → prim ('&'? prim)*             -- highest; '&' implicit
prim          → '!'? (isotope | element_symbol | hydrogen | wildcard | chirality
                     | charge | atomic_number | property_primitive
                     | recursive | atom_class)
recursive     → '$(' chain ')'
```

Four things this fixes that a direct spec transcription gets wrong, all found by
running `tree-sitter generate`:

- **No rule may match the empty string.** `side → chain?` is a hard generate
  error, not a warning. Express the empty agent side of `>>` at the use site.
- **`chain`'s tail is `optional`, not `repeat`.** `branched_atom (bond? chain)*`
  is ambiguous — `CCC` splits two ways — and fails to generate. A single optional
  recursive tail is unambiguous and is the form that generated cleanly in a spike.
- **`conflicts: $ => [[$.branched_atom]]` is required.** After an atom, `-` may
  begin a ring closure (`C-1CCCCC-1`, valid) or a chain bond (`C-C`). LR(1)
  cannot decide.
- **`charge`, `chirality`, `atom_class`, `atomic_number` and
  `property_primitive` must each be `token(...)`.** Otherwise implicit `&`
  between primitives makes `[--]` (charge −2 vs two negations), `[#6]` (`#6` vs
  `#` & `6`) and `[C:1]` ambiguous.

**Logical levels are named for their operator.** SMARTS precedence is `&` > `,` >
`;`, so the outside-in nesting is low-AND(`;`) → OR(`,`) → high-AND(`&`). Naming
the `;`-operand "or_expr" and the `&`-repeat "and_high" keeps each rule name
describing its own body; every query and corpus file inherits those names.

**`!` is a prefix on any bond primitive**, not a special `!@` token: `C!-C`,
`C!=C`, `C!#C` all parse, and `C!-!=C` reads as `C!-&!=C`.

**`recursive` nests a `chain`, not a `structure`.** `[$(C>>C)]` is rejected;
`[$(C.C)]`, `[$([$(C)])]` and `[$([C:1])]` all parse.

**`atomic_number` is its own node**, separate from `property_primitive`: `[#7]`
is nitrogen, so the element-colouring layer has to reach it.

**`chirality` covers the extended classes**, not just `@`/`@@`:
`'@' ('@' | ('TH'|'AL'|'TB'|'OH') digit+)?`, all four measured as parsing.

**`hydrogen` is not in the element table.** With `H` present as a string literal
it beats a `/H[0-9]*/` regex at equal length, and `[nH]` — among the most common
brackets in the language — parses as two element symbols instead of `n` plus an
H-count. Give hydrogen its own token with explicit lexer precedence
(`token(prec(1, /H[0-9]*/))`) and decide `[H]`-alone deliberately.

**`property_primitive` is a probed set, not a guess.** `D d h R r v X x z Z k ^`
all parse (`[c;h1]` → `[c&h1]`, `[C^3]` → `[C&^3]`); `h` and `^` are in everyday
use, and omitting one means ERROR nodes on valid input. Enumerate from probing,
not from the spec.

**Generate the symbol table from the parsers, not the periodic table.**
`GetPeriodicTable()` yields 119 entries, and Z=0's symbol is `*` — which would
land in the element list and collide with the wildcard token. The accepted SMARTS
symbol set is also a strict subset of the elements (`Mc`, `Ts`, `Og`). So probe
`MolFromSmarts("[" + sym + "]")` and `MolFromSmiles(...)` per symbol and emit the
four disjoint sets the parsers actually accept: bare aliphatic, bare aromatic,
bracket element, bracket aromatic.

**CXSMILES starts shallow**: recognise the `|…|` layer, its `$…$` label block,
and `key:value` fields generically. The field vocabulary (`atomProp`, `Sg`, `m`,
`^`, `w`, `c`/`t`, `ha`, coordinates) is a sublanguage of its own and can deepen
later without reshaping anything above it.

## Layout and wiring

`modules/tree-sitter-smarts/`, mirroring `tree-sitter-miller`: `grammar.js`,
`symbols.js`, `symbols.js.py`, `package.json`, `tree-sitter.json`,
`queries/highlights.scm`, `test/corpus/`. Commit `src/` — nvim-treesitter builds
from it and only regenerates when `install_info.generate` is set, so a second
machine needs no `tree-sitter generate`. Run `generate` with `tree-sitter.json`
present so the committed `parser.c` carries `LANGUAGE_VERSION 15`, matching
`vim.treesitter.language_version`. `*.so` is already gitignored; the build lands
`parser.so` in the module directory.

Queries flat at `queries/highlights.scm`, with `tree-sitter.json`'s `highlights`
field pointing there and `injection-regex` set to the grammar name.

- `plugin/ftdetect.lua` — `smi`, `smiles`, `cxsmiles`, `rsmi` → `smiles`;
  `smarts`, `sma`, `smirks` → `smarts`.
- `lua/plugins/treesitter.lua` — `parsers.smarts` with `install_info = { path = …
  .. '/modules/tree-sitter-smarts', queries = 'queries' }` in the `User TSUpdate`
  callback, and `"smarts"` in the `install` list. That whole block is gated on
  `vim.fn.executable("tree-sitter")`, so a machine without the CLI gets no parser
  and no highlighting.
- `plugin/treesitter.lua` — `language.register("smarts", "smiles")`, alongside the
  existing registrations. This also makes ```` ```smiles ```` fences resolve,
  since `language.get_lang` consults the same alias table.
- `ftplugin/smiles.lua` and `ftplugin/smarts.lua` — `commentstring = "#%s"`,
  matching this repo's convention of no space. One structure per line means no
  indenting and no `indents.scm`. While nothing else diverges, the second can
  `runtime! ftplugin/smiles.lua`.

## Highlights

Structural captures: bonds `@operator`, branch parens `@punctuation.bracket`,
brackets `@punctuation.special`, ring closures `@number`, isotope `@number`,
charge `@operator`, atom class `@label`, chirality `@attribute`, SMARTS logical
`,;&!` `@keyword.operator`, `property_primitive` `@function.builtin`,
`atomic_number` `@number`, `$(` `@punctuation.special`, the name field `@string`,
the header `@title`, comments `@comment`.

Elements stay one node type per context — `organic_symbol`, `aromatic_symbol`,
`element_symbol` — with no per-element distinction in the tree. Which element a
node is, is its text, and text matching belongs in a query. Base captures give
aliphatic one group and aromatic another.

### Element colouring layer

A second query layer matches element identity by text and assigns a private
group. It has to span all three element node types, plus `atomic_number` for the
`[#7]` spelling:

```scheme
; extends
([(organic_symbol) (aromatic_symbol) (element_symbol)] @chem.nitrogen
 (#any-of? @chem.nitrogen "N" "n")
 (#set! priority 101))

((atomic_number) @chem.nitrogen
 (#eq? @chem.nitrogen "#7")
 (#set! priority 101))
```

- `#any-of?` tests the captured node's own text — `:help
  treesitter-predicate-any-of?`.
- Explicit `priority 101` over treesitter's default 100 (`:help
  treesitter-highlight-priority`), rather than relying on pattern order within the
  file.
- `@chem.*` is a private namespace, so no colorscheme paints it by accident.
  Define the groups in `hi.onColorScheme`.
- Lives at `after/queries/smarts/highlights.scm` in the nvim config, not in the
  grammar module: the module stays free of this config's colour opinions, and the
  whole layer toggles by renaming one file while comparing palettes.

First palette to try is CPK — N blue, O red, S yellow, P orange, halogens green,
C the default foreground — judged on rendered examples.

## Tests

`test/corpus/`, run with `tree-sitter test`. The cases that carry the bugs:

- Element alphabets: `Sn` as two atoms vs `[Sn]` as one, `Sc1ccccc1`, `[Ra]`, and
  `[Nh]` as the documented divergence.
- Brackets: `[nH]` and `[CH]` (H with no digit), `[--]` and `[++]`, `[#6]`,
  `[C:1]`, `[c;h1]`, `[C^3]`, nested `$( )`.
- Bonds and closures: `C$C`, `C!-!=C`, `%(100)`, bond-prefixed closure
  `C-1CCCCC-1`.
- Reactions: `>>` with an empty agent side, and a three-side `C>C>C`.
- Record layer: comment line, blank line, header line, trailing whitespace after
  a structure, CRLF, a spaced name, a name beginning with `|`, and a CX layer on
  both a molecule and a reaction.

`tests/plenary/smarts_spec.lua` for the nvim side: extension → filetype for each
extension (`.smi` in particular, against `mib`), filetype `smiles` resolving to
the `smarts` parser, and a parse smoke test asserting no ERROR node over a small
file. Per `tests/README.md`, call `get_parser(buf):parse()` before inspecting.

## Follow-ups

Each is independent of the others and needs the grammar stable first.

1. **Markdown and typst fenced blocks** — free once the parser exists, since the
   fence info string supplies `injection.language`. Confirm with a fixture.
2. **Python injection at rdkit call sites** — a query matching
   `MolFromSmiles`/`MolFromSmarts`/`ReactionFromSmarts` and capturing the string
   argument. Static and precise, unlike the buffer-wide import sniffing in
   `ftplugin/python.lua`.
