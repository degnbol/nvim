# DRAFT PLAN: argparse attribute diagnostics

Catch `args.<typo>` — access of a namespace attribute that no `add_argument`
ever created. basedpyright cannot: typeshed gives `argparse.Namespace` a
`__getattr__(self, name: str) -> Any` (`argparse.pyi:493`), so every attribute
access type-checks. This is deliberate in typeshed and no pyright setting changes
it — pyright has no plugin hook to teach it argparse's runtime dest derivation.

So the diagnostic has to come from outside the LSP: a small linter that reads the
buffer, derives the set of dests the parser will actually create, and flags
`args.X` accesses not in that set.

**This is the typo class only** (`args.term` when the dest is `args.query`). The
separate "dest silently comes from the *first* long option, not the longest"
footgun is **out of scope** — it's not an error, and flagging every multi-long
`add_argument` is noise the user explicitly rejected.

## Scope and bail-out (false-positive control)

The check is sound only for a **flat, single-parser, single-`parse_args`** buffer.
These constructs move dests out of view or make them conditional, so their
presence must **suppress all diagnostics for the buffer** (emit nothing) rather
than risk false positives:

- `add_subparsers()` — each subparser owns a separate dest namespace; dests are
  conditional on the chosen subcommand.
- `parents=[...]` — dests defined on another parser object, often another module.
- `set_defaults(**kw)` — injects dests with no `add_argument` call.
- `getattr(args, name)` / `**vars(args)` / `argparse.Namespace` passed around —
  dynamic access the linter can't resolve.

Bail-out beats false positives: a linter that cries wolf on colleagues'
subparser-based CLIs gets disabled. Silence on the hard cases, accuracy on the
scripts where the user actually makes the typo (flat uv single-file scripts).

## Why Python `ast`, not `vim.treesitter` (even though it's native)

`vim.treesitter` is in-process, zero-dependency, and already has the buffer's
parse tree for free — the obvious default. It is nonetheless the wrong tool here,
for reasons specific to argparse:

1. **Dest derivation is string logic over argument *values*, not tree structure.**
   argparse computes a dest by: taking the option strings, choosing the first
   `--long` (else the first `-x`), stripping leading `-`, and replacing `-`→`_` —
   unless `dest=` overrides it. treesitter locates the `add_argument` node and its
   string arguments, but the derivation itself is imperative string munging you'd
   reimplement in Lua regardless. treesitter buys nothing for the actual hard part.

2. **`ast` can reach ground truth; treesitter categorically cannot.** Because an
   `ast` linter runs *in Python*, it has an escalation path treesitter has no
   analogue for: build the parser and read `parser._actions[*].dest` — argparse's
   own answer, with zero reimplementation of the first-long-option rule. (Guarded:
   only viable if parser construction is side-effect-free; module-level
   `parse_args()`/`sys.exit` means fall back to static `ast`. Still, the ceiling
   exists only on the Python side.) treesitter can only ever pattern-match source
   text; it can never *ask argparse*.

3. **Keyword-argument inspection is first-class in `ast`.** `dest=`, the
   `parents=`/`set_defaults`/`add_subparsers` bail-out triggers — all are
   `ast.keyword` nodes with `.arg`/`.value`, trivially read. In treesitter you walk
   `keyword_argument` nodes matching identifier text by hand.

4. **Reuse and single-source-of-truth — the decisive one.** This config already
   shares Python type infrastructure between neovim and the Claude lint hook
   (`lsp_ext/` stubs consumed by both `lsp/basedpyright.lua` and
   `lib/lint-tier.sh`). A Python `ast` script is reusable the same way: nvim
   diagnostics, the lint hook, pre-commit, CI. A `vim.treesitter` query is locked
   inside nvim and can never back the shell lint hook. Global CLAUDE.md:
   "Deterministic over probabilistic — a parser over regex"; here the most
   authoritative parser *is* argparse itself, reachable only from Python.

Honest counter-argument: for a keystroke-frequency in-editor check, native
treesitter avoids the subprocess spawn. It doesn't apply — nvim-lint (and a
diagnostic-on-save autocmd) spawns a subprocess by design, and the linter should
run on save/`BufWritePost`, not per keystroke. The subprocess cost is paid once
per save and buys reusability + the ground-truth ceiling.

## The linter script

`lsp_ext/argparse_dests.py` (see §"Location" below), uv shebang, stdlib-only
(`ast`, `sys`, `json`):

1. Read source from stdin (or `$1`). `tree = ast.parse(src)`.
2. Walk for bail-out triggers: any `Call` to `.add_subparsers`, any `add_argument`
   with `parents=` (actually `parents=` is on `ArgumentParser(...)`), any
   `.set_defaults(...)`, any `getattr(args, ...)`. If found → print `[]`, exit.
3. Collect dests: for each `.add_argument(...)` call, derive dest from `dest=` kwarg
   if present, else from the option-string positional args via the rule in §2.1.
   Positional-argument `add_argument("name")` → dest is the name verbatim.
4. Identify the namespace variable: LHS of `x = <parser>.parse_args(...)`
   (and `.parse_known_args()` → tuple, first element). Usually `args`.
5. Walk `Attribute` accesses `<nsvar>.X` in `Load` context; for each `X` not in
   the dest set, emit `{lnum, col, end_col, message}` (0-indexed for
   `vim.diagnostic`). Skip `Store`/`Del` context and dunder/known-real attrs
   (`Namespace` provides none worth allowing beyond the derived set).
6. Print JSON array to stdout.

Standard-library dests to seed the known set regardless: none — argparse adds no
implicit dests except `help`, which is an action with no attribute. So the derived
set is complete for the flat case.

## Wiring into neovim

No lint plugin is installed (`nvim-lint` absent; `lua/plugins/` has no lint spec).
Two options:

- **A — plain autocmd + `vim.diagnostic` (no new dependency).** `BufWritePost`
  (and `BufReadPost`) on `python` filetype → run the script via `vim.system`,
  `vim.json.decode` the output, `vim.diagnostic.set(ns, bufnr, items)` under a
  dedicated namespace `vim.api.nvim_create_namespace("argparse_dests")`. Matches
  the config's existing custom-diagnostic pattern (`lua/tex/textcolor.lua`,
  `lua/plugins/lsp.lua:87`). Lives in `ftplugin/python.lua` or a new
  `lua/lint/argparse.lua` required from there. **Recommended** — one script, one
  namespace, zero plugin surface.
- **B — adopt `mfussenegger/nvim-lint`.** Register the script as a custom linter
  with a parser. Heavier (new plugin, lz.n spec, pack_specs entry) and only worth
  it if other ad-hoc linters are wanted later. Defer unless that need is real.

Diagnostic severity: `WARN` (it's a likely typo, but the bail-out means a
surviving false positive shouldn't hard-error a lint gate).

## Location: `lsp_ext/`

The single copy lives at `lsp_ext/argparse_dests.py` — the established
shared-between-nvim-and-lint-hook Python-tooling directory, which makes the reuse
argument (§3.4) real: both `ftplugin/python.lua` and `lib/lint-tier.sh` can reach it.

Two clarifications so the fit isn't overstated:

- **It's an executable, not passive type data.** Everything else in `lsp_ext/`
  (`python_stubs/`, `extraPaths/`) is consumed *by pyright* to extend its inputs.
  This script instead *produces* diagnostics pyright structurally can't. The
  precedent for an executable in `lsp_ext/` is `r_lsp_dots.R` — this is a sibling
  to it, not to the stub dirs.
- **No automatic wiring.** The "drop it in and both consumers glob it" property is
  specific to `extraPaths/` for import resolution. This script must be called
  explicitly from both `ftplugin/python.lua` (the nvim diagnostic autocmd) and
  `lib/lint-tier.sh` (if the lint hook is to surface it too). Co-location gives the
  shared home; the calls are still hand-wired.

## Open questions

1. **Runtime introspection (§3.2) — build it now or later?** Static `ast` covers
   the typo class fully for flat parsers. The `parser._actions` escalation only
   matters for exotic dest derivation (custom `action=` classes changing dest).
   Recommend static-only v1; note the escalation path but don't build it.
2. **`parse_known_args` and reassigned namespaces** (`args = parser.parse_args();
   args = fixup(args)`) — v1 can track only the direct `parse_args` LHS and bail if
   the namespace var is reassigned from a non-parse_args expression.
