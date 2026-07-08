# PLAN — `grd` goto-def for glossy `@term` refs

Make `grd` in typst buffers jump to a glossary entry's *definition* (the dict
pair in the source `.typ`) for glossary terms — `@ts`, `@ts:pl`, `@ts:short`, …
— while leaving `@fig:`/`@sec:`/`@eq:` on the (correct) LSP path.

Mechanism menu: `~/.claude/skills/neovim/references/keymaps.md`.

## Why the LSP alone is wrong here

Verified via probes (see typst-ref.md):

- `@fig:f` / `@sec:h` / `@eq:e` → LSP resolves to the `<label>`. **Correct. Keep.**
- `@ts` (glossy term) → LSP lands on `glossy/gloss.typ:624`, the package's
  runtime `label(key)` placeholder. **Useless.**
- `@ts:pl` (modified) → LSP resolves to nothing (`<ts:pl>` label never created).

So precedence is per-input-class: the custom handler must **win** for glossary
keys and **defer** for everything else. "LSP first, custom fallback" is wrong —
the LSP returns a valid-but-useless hit for `@ts`, which would shadow us.

## What a glossary entry looks like

Entries are Typst dict pairs assigned to a var and merged in the preamble via
degnlib's `merge-glossary` → `make-glossy-init`, e.g.
`~/Documents/xyme/typst/reactions/glossary.typ`:

```typst
#let glossary-local = (
  ts: ( short: "TS", long: "transition state", description: [...] ),
  api: ( short: "API", ... ),
)
```

The `@`-ref key is exactly the dict key (`ts`, `api`). Definition site = the
line `  ts: (`. Special-char keys are quoted (`"1-2-shift": (`, `"π-bond": (`).

## Layout / where the sources are

The document declares its glossary sources by `#import`, all relative paths:

```typst
main.typ:      #import "preamble.typ": *
preamble.typ:  #import "./glossary_base.typ": glossary-base
               #import "./glossary.typ": glossary-local
               #let glossary-data = merge-glossary(glossary-base, glossary-local)
```

`glossary_base.typ` is a per-doc symlink to the shared base one dir up. **The
import statements are the specification** — filenames (`glossary.typ`,
`glossary_base.typ`, could be `glossary_local.typ`) are convention. So discovery
follows the import graph, never guesses filenames or walks the filesystem.

## Architecture

Pure logic in a require-able module `lua/typst_glossary.lua`; `ftplugin/typst.lua`
only wires the buffer-local keymap. This is what makes the logic testable in
isolation (see Check).

Module API:

- `ref_key(text)` — ref node text → key.
- `imports(content, dir)` — treesitter → resolved relative import paths.
- `entries(content)` — treesitter → `{ key = row }` for glossary dict pairs.
- `resolve(bufnr)` — BFS the imports, match the key → qf items (or `nil`).

## Resolution — treesitter throughout

Every step uses the typst treesitter grammar (parser already installed; node
types verified against the real files). No Lua line-patterns, no `rg`, no
`vim.fs.find`.

1. **Ref key.** `vim.treesitter.get_node` under cursor. If not a `ref` (or ref
   ancestor) → defer to LSP. Else key = ref node text minus leading `@`, before
   the first `:` (`text:gsub("^@",""):match("^[^:]+")`). This yields `ts` for
   `@ts`/`@ts:pl`/`@ts:short` and `fig` for `@fig:f`. A `gsub` on the node's own
   text, not a buffer scan. Robust across the whole ref span incl.
   hyphenated/unicode keys (`@wagner-meerwein`, `@π-bond`).

2. **Discovery — BFS the relative import graph.** From the buffer, query imports
   and recurse (visited set), resolving relative paths against each importing
   file's dir:

   ```scheme
   (import (string) @path)
   ```

   Keep paths starting with `.` (relative); skip `@preview/…` and `@local/…`
   package imports. Reads only files the document actually imports, so sibling
   documents are unreachable **by construction** and package deps are ignored.
   The graph is a handful of tiny files — parse per `grd` press, no caching.

3. **Entry match.** For each reachable file, `get_string_parser(readfile)` and
   query:

   ```scheme
   (tagged [(ident) (string)] @key (group)) @entry
   ```

   The `(group)` value constraint matches entry-level pairs `key: (…)` (bare or
   quoted key) and excludes string-valued fields (`short: "TS"`). Filter
   captures where `@key` text (quotes stripped) == the ref key; take the node's
   row. Build qf items `{ filename, lnum, col, text }`.
   - **Hits → `map.qf_mini({ items = items })`** (jumps if 1, quickfix if
     several; preserves jumplist via `util.jump`).
   - **No hits → defer to LSP** (key isn't a glossary entry, e.g. `fig`/`sec`).

   Known collision surface: the query also matches any `key: (group)` in a
   reachable non-glossary file (e.g. a config dict `margin: (…)` in
   `preamble.typ`). Requires a config field named exactly like a glossary term;
   harmless (2-item quickfix). Scope to `merge-glossary`'s arg imports only if it
   ever bites.

4. **LSP defer.** Regular buffer-local map, so on miss just call the LSP
   definition directly via a shared helper. The global `grd`
   (`plugin/keymaps.lua:384`) and this map both want the same
   `vim.lsp.buf.definition(filter_lsp_items(…))` body — extract it to
   `utils/keymap.lua` as `M.lsp_definition()` and call it from both. No
   `maparg` capture-and-chain: that mechanism is for wrapping mappings you don't
   own; we own the global, so we refactor instead. (Note `ftplugin/typst.lua`'s
   `typst_lsp` reference is dead — the active server is **tinymist**,
   `lsp/tinymist.lua`; removed below.)

## Install location

Buffer-local `grd` in `ftplugin/typst.lua`
(`map.n('grd', fn, "...", { buffer = true })`), calling `require"typst_glossary".resolve`
and `map.lsp_definition()` on miss. Buffer-local shadows the global — no
fallthrough, so the miss branch calls the helper explicitly.

## Dead code removed in the same change

Delete the `<LocalLeader>t` "toggle update on type" map in `ftplugin/typst.lua`
(lines 13–27): it calls `require"lspconfig".typst_lsp`, defunct since the
tinymist migration, via the outdated `lspconfig.setup` + `LspStop`/`LspStart`
pattern. On-type PDF export is a *setting* — if ever wanted, set
`exportPdf = "onType"` in `lsp/tinymist.lua`, not a runtime toggle.

## Check to leave behind

Plenary spec `tests/plenary/typst_glossary_spec.lua`, extending the existing
`describe` style (`tests/README.md`), three fixture-free pure tests:

1. `ref_key`: `@ts`/`@ts:pl`/`@ts:short` → `ts`; `@fig:f` → `fig`;
   `@π-bond` → `π-bond`.
2. `entries` on an inline dict string: finds `ts` and dequoted `1-2-shift` with
   correct rows; does **not** find `short` (string-valued field).
3. `imports` on an inline preamble string: keeps `./glossary.typ`, drops
   `@preview/…` and `@local/…`.

Skip an on-disk BFS integration test: sibling docs are unreachable by
construction (only explicit relative imports are followed), so it would defend a
property the design already guarantees.

## Open items to confirm during build

- A buffer with no glossary up its import graph (scratch/nameless, or a stray
  `.typ`) — BFS yields nothing → defer to LSP on the empty result.
- Editing a glossary file directly: it isn't imported by itself, but `grd` there
  targets `@refs` in prose, which don't occur in the dict file — non-case.
