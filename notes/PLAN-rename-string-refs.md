# DRAFT PLAN: extend LSP rename to structural string references (Tier 1)

## Goal

When renaming a symbol, also update string literals that *structurally* refer to
it but which the LSP treats as opaque text. Motivating case (Tap / argparse):

```python
class Args(Tap):
    infiles: list[str]
    out: str = "/dev/stdout"

    def configure(self) -> None:
        self.add_argument("infiles", nargs="+")   # "infiles" missed by rename
        self.add_argument("-o", "--out")          # "--out" missed by rename
```

Renaming `infiles` → the positional string `"infiles"` must change too; renaming
`out` → the `--out` option string (and a `dest="out"` if present) must change.
`-o` (short) stays.

**Only structural references.** Not "any string equal to the name" — match the
call *shape* (`add_argument(...)` whose argparse-derived dest equals the old
name). A coincidental `"infiles"` elsewhere is never touched. This keeps Tier 1
deterministic and confirmation-free. Prose/comments are the separate, non-editing
Tier 2 — see [PLAN-rename-illuminate-scope.md](PLAN-rename-illuminate-scope.md).

## Why the server can't do it

basedpyright has no static link from the `infiles` attribute to `"infiles"` —
Tap builds the dest at runtime, and pyright has no plugin hook (same root cause
documented in [PLAN-argparse-diagnostics.md](PLAN-argparse-diagnostics.md)). Must
be client-side.

## The hook: a per-client rename handler wrapper

The single choke point where *any* rename applies its `WorkspaceEdit`:

- live-rename's `submit()` calls
  `ctx.client.handlers['textDocument/rename'] or vim.lsp.handlers['textDocument/rename']`
  with `(err, result, context, config)` — `live-rename.lua:738-740`. `result` is
  the server's `WorkspaceEdit`.
- Stock `vim.lsp.buf.rename` resolves the handler identically (`buf.lua:765`).

So install, at `LspAttach` (in `lua/autocmds/lsp.lua`, gated to the relevant
filetypes/clients), a wrapper:

```lua
local default = vim.lsp.handlers['textDocument/rename']
local wrapped = function(err, result, ctx, config)
    if not err and result then
        augment_workspace_edit(result, ctx)  -- mutate in place, add TextEdits
    end
    return default(err, result, ctx, config)
end
wrapped._tier1 = true  -- sentinel, see below
client.handlers['textDocument/rename'] = wrapped
```

This augments *before* the one apply, so everything lands in a single undo block
and there is no double-apply / offset-shift bookkeeping. It also decouples Tier 1
from the `grn` keymap: the wrap is agnostic to which rename UI is bound (live-rename
or stock `vim.lsp.buf.rename`), and any direct `textDocument/rename` caller gets it.
(Not code actions — those apply their `WorkspaceEdit` straight through
`apply_workspace_edit`, `buf.lua:1259-1260`, never touching the rename handler. A
rename delivered as a code action is out of scope.)

**Install once per client, not per attach.** `LspAttach` fires once per
(client, buffer) pair, but `client.handlers` is a single shared object on the
client. Re-running the wrap on every attach would wrap the already-wrapped
handler, so `augment_workspace_edit` runs 2×, 3×… emitting duplicate TextEdits.
Guard on the sentinel (skip if `client.handlers['textDocument/rename']._tier1`)
or track installed `client_id`s.

### Deriving old/new name inside the handler

- **new name**: `ctx.params.newName` (rename params carry it), or the `newText`
  of any edit.
- **old name**: read the current buffer text at one edit's range *before* the
  default handler applies it (every edit replaces old→new, so any range yields
  the old text). Do this per-file for the files the edit touches.
- **Load target buffers first.** For a multi-file WorkspaceEdit, reading the old
  text and running the matcher's treesitter parse both require the target buffer
  to be *loaded*. `vim.uri_to_bufnr` creates the buffer but does not load it, and
  we run *before* `apply_workspace_edit` (which would load it). `fn.bufload(bufnr)`
  each target URI before reading/parsing.

### The trust boundary

Only augment files/ranges the server's `WorkspaceEdit` already includes. The
server decided those files are part of this symbol's rename; we add string edits
*within* them. Never scan the whole workspace independently.

## The matcher registry

A table of matchers keyed by a predicate (filetype + AST shape). Each matcher,
given a file + old name + new name, returns extra `TextEdit`s. Ship one:

**argparse/Tap `add_argument`** (treesitter over the file's tree):
1. Find `call` nodes whose function is `add_argument` (attribute or bare).
2. Derive the dest from the string arguments using the **same logic** as
   [PLAN-argparse-diagnostics.md](PLAN-argparse-diagnostics.md): first `--long`
   (else first `-x`) with leading dashes stripped and `-`→`_`, unless a `dest=`
   kwarg overrides. For a positional (no leading dash), dest = the name itself.
3. If derived dest == old name, emit a `TextEdit` over the *specific* string
   token that produced the dest — the flag edit and the `dest=` edit are
   **mutually exclusive**, never both:
   - **Positional** (no leading dash): rewrite the name string.
   - **Dest derived from a flag** (no `dest=` kwarg): rewrite the `--long` option
     it came from. dest-derivation maps flag `-`→`_`, so invert it when writing
     the flag: the new flag is the new name with `_`→`-` (`out`→`--out`,
     `out_file`→`--out-file`). This round-trips through argparse's own
     normalization and matches CLI convention. This edit *is* required — for
     Tap/argparse the flag and the attribute are coupled, so renaming the
     attribute without the flag breaks the binding.
   - **Explicit `dest=`**: rewrite only the `dest=` value. The option strings are
     independent user-facing CLI surface here and must **not** be touched —
     `add_argument("--out", dest="out")` renamed `out`→`output` edits `dest`
     only, leaving `--out` alone.

The matcher fires on the structural gate alone — no restriction to `Tap`
subclasses or files importing argparse/tap. A stray `x.add_argument("infiles")`
with a real matching dest is correctly a reference, and an import check adds a
fragile heuristic that buys nothing.

Future matchers (not now): `getattr/setattr/hasattr` literals, `pydantic
Field(alias=)`, dataclass field strings. The registry shape is the deliverable;
one matcher proves it.

## Non-obvious requirements

- **Offset encoding.** Added `TextEdit` ranges must be in the client's
  `offset_encoding` (utf-16 by default), same as the server's edits. Compute
  columns accordingly (`vim.lsp.util` / `vim.str_utfindex`), not raw byte
  columns.
- **Non-overlap.** Server edits cover identifier occurrences; ours cover string
  literals — disjoint by construction. `apply_workspace_edit` sorts and applies;
  just don't emit an edit onto a range the server already has.
- **Both `changes` and `documentChanges` shapes.** The `WorkspaceEdit` may use
  either (live-rename's own ref handler already branches on both,
  `live-rename.lua:208-222`). Augment whichever is present.
- **Structural gate, not text gate.** Matching is call-shape + derived-dest
  equality. Do not fall back to "string literal equals old name" — that reopens
  the false-positive door Tier 1 exists to avoid.

## Shared dest-derivation with the argparse diagnostic

[PLAN-argparse-diagnostics.md](PLAN-argparse-diagnostics.md) implements the same
first-long-else-short / strip-dashes / `dest=`-override rule as a Python `ast`
script for the shell lint hook. Tier 1 needs the same rule in-process (Lua over
treesitter) — a subprocess in the interactive `grn` hot path is unacceptable, and
Tier 1 needs source *token ranges* the diagnostic's runtime approach can't give.
So: **keep two implementations, single-source the *truth* via a shared fixture
set** both must agree on in tests (a table of `add_argument(...)` snippets →
expected dest). This is the deterministic guarantee without dragging a subprocess
into rename. Not an open question — decided.

## Module layout

A new `lua/lsp_rename/` module (handler wrapper + registry), required from
`lua/autocmds/lsp.lua`. Matchers are pure `(src, old, new) -> TextEdit[]`,
unit-testable without a live server.

Most of the module is rename/argparse-specific and stays. The generic LSP
primitives it needs already exist in `vim.lsp.util`: offset-encoding column
conversion is `character_offset(buf, row, col, offset_encoding)`, range params
`make_given_range_params`. The one primitive the runtime lacks is normalizing a
`WorkspaceEdit` to a shape-agnostic view: given a `WorkspaceEdit` (either
`changes` or `documentChanges`), return the per-file `{ uri, edits }` pairs where
each `edits` is the live `TextEdit[]` array. That is its single job — the caller
appends its own `TextEdit`s into the returned array with `table.insert`. (The
apply side handles both shapes, but there's no public normalizer — live-rename
reinvents it at `live-rename.lua:208-222`.) General LSP plumbing with no rename
coupling, so it lives in a new `lua/utils/lsp.lua` — placed there for its
generality, not gated on a second consumer.
