# Docstring reconcile — Plan

A docstring **reconciler**: keep an existing docstring in sync with its function
signature — add missing parameters, drop into signature order, preserve
hand-written descriptions. Cross-language by design; Python/Google first.

This is the gap no existing tool fills. codedocs / neogen / vim-doge /
vim-python-docstring / docstring.nvim / neodoc.nvim are all *generate-only*
(surveyed 2026-08-07; neogen's sync request #65 is open/help-wanted, vim-doge's
#18 is wontfix). pyment (CLI) is the only tool that preserves descriptions on
regenerate, but it neither reorders nor prunes.

Clean-room implementation — codedocs is GPL-3, so we take *ideas* (style-as-data,
treesitter extraction, golden-file fixtures) but no code, leaving us free to
license as we choose.

## Where this lives — core `lua/`, not a module

The docstring parse/render/model layer is **general infrastructure**, not one
feature's private domain: it is consumed by reconcile *and* by an `lsp_rename`
docstring matcher, and could serve future consumers (a hover formatter, a
docstring diagnostic). By name-for-scope it is a *docstring* library, not
"doccer's" library, so it belongs in core `lua/docstring/`.

Once the library is core, the reconciler is just **another consumer** of it,
exactly like the rename matcher. Isolating that one consumer in a
`modules/*.nvim` dev plugin would only manufacture a cross-boundary dependency
(module → core `lua/docstring/`) that makes the "module" non-independent for no
gain — the code is never published and needs no load isolation. So the whole
thing is core `lua/`; `modules/` stays reserved for genuinely separate repos.

```
lua/docstring/
  model.lua            -- DocModel type
  style/<style>.lua    -- parse / render / locate_doc_entries  (google, numpy, rest)
  lang/<language>.lua  -- treesitter: locate_target / signature_params / locate_docstring
  reconcile.lua        -- the reconcile algorithm (a consumer)
  command.lua          -- :Doccer-style command (a consumer)
lua/lsp_rename/docstring.lua   -- rename matcher (the other consumer)
```

Everything points at `lua/docstring/`; the consumers never depend on each other.

## Scope

**In scope**

- Add parameter (at end or middle) → insert entry in the correct position.
- Reorder parameters → reorder doc entries to signature order.
- Preserve hand-written descriptions on unchanged params.
- Generate from scratch when there is no existing docstring.
- Idempotent: `reconcile ∘ reconcile = reconcile` (see below).

**Out of scope**

- **Rename.** Renaming a parameter belongs in `lsp_rename` (see [Docstring
  rename](#docstring-rename-belongs-in-lsp_rename)). A name in the docstring but
  absent from the signature is ambiguous (renamed vs removed) and can't be
  resolved from the signature alone; such "orphan" entries are *surfaced, never
  silently dropped*.
- Rewriting or improving prose descriptions.
- Class-attribute docstrings (later milestone).

## Why from scratch, not fork

codedocs models the **signature** language-agnostically (`item_extractor` →
`{name, type}` via per-language treesitter queries) and the **output** as style
`blocks`/`layout`. It has **no reverse direction** — nothing parses an existing
docstring back into structure. Its render path is `vim.snippet.expand`
(interactive generate with tabstops), architecturally opposed to reconcile,
which must compute final text and *replace a range* without re-prompting for
preserved content. So the reusable surface is small and the hard half is new
regardless.

## Architecture

Three transforms, all language/style-parameterised:

```
signature ──signature_params──▶ [SigParam]   (forward: treesitter, per language)
docstring ──parse────────────▶ DocModel      (reverse: prose→structure, per style)  ← NEW
DocModel  ──render───────────▶ text          (forward: structure→prose, per style)
```

Reconcile lives between them:

```
reconcile(SigParam[], DocModel) ▶ DocModel'
```

### DocModel — the language-agnostic object

A parsed docstring as ordered sections, preserving anything we don't understand
so we never lose data:

- `summary` — free text (untouched).
- `params` — ordered list of `{ name, type?, description }` (description may be
  multiline).
- `returns` — `{ type?, description }` (nameless — single slot).
- `raises` — list of `{ type, description }`.
- `other` — any section we don't model, kept verbatim in place.

### Pipeline (`:Docstring` on cursor)

1. Locate the target node under the cursor (function/method) — treesitter.
2. Extract ordered `SigParam[]` from the signature — treesitter (filtering rules
   below).
3. Find the existing docstring node + its buffer range — treesitter.
4. No docstring → render a fresh DocModel from the signature.
5. Docstring exists → `parse` → `reconcile` → `render` → replace the range.

### Signature extraction rules

`signature_params` is not a raw child list; the Google convention and Python
grammar require filtering, decided up front (milestone 1):

- **`self` / `cls`** — omitted from `Args:`, so dropped.
- **`*args` / `**kwargs`** — kept, documented under those names (with the stars)
  as ordinary entries.
- **`/` and `*` markers** (positional-only / keyword-only separators) — not
  params; skipped.
- **defaults** — do not affect the entry name; noted only if a style renders
  them.
- **empty signature** → fresh-gen produces a docstring with no `Args:` section.

### Reconcile algorithm (params)

- Order comes from the signature.
- For each `SigParam`, reuse the matching doc entry by name and keep its
  description. **Type**: use the signature annotation when present (code is the
  source of truth); when the signature param is *unannotated*, **preserve the
  docstring's existing type** — never blank a hand-written type just because the
  code lacks an annotation.
- Doc entries with no matching signature name = **orphans** → held in an
  `unmatched` bucket, rendered under a visible marker rather than dropped. The
  marker format must itself **round-trip**: a second reconcile must re-parse it
  back into `unmatched`, not into `other` or a duplicate (see Testing).
- `returns` — single nameless slot: preserve the description, refresh the type
  from the return annotation.
- `raises` — **verbatim passthrough**. Nothing in a signature enumerates raised
  exceptions (they are body-level `raise` statements), so there is no source to
  refresh from. (Phase 7.)

**Idempotency** is `reconcile ∘ reconcile = reconcile`, not text-identity: the
first run may normalise a semantically-correct but hand-formatted docstring
(type placement, whitespace). The guarantee — and the test — is that the *second*
run is a no-op.

### Buffer edit

- `nvim_buf_set_text` over the docstring node's range with the rendered text.
  Multiline rendered text is indented **per line** to the function body indent,
  not inserted as one blob. **Not** snippet expansion.
- Phase 2 nicety: place snippet tabstops on *new blank* entries only, so the
  user can fill them, while preserved text stays inert.

## Docstring rename belongs in lsp_rename

Renaming a parameter is a *symbol* operation: the LSP renames the identifier and
its body uses, but treats the docstring's `Args: name: ...` line as opaque text
and leaves it stale. That is precisely the class of problem this config's
`lsp_rename` already solves for argparse/Tap option strings — structural string
references the server ignores, augmented into the server's `WorkspaceEdit`
*before* it applies so the string edit lands in the same undo block.

So docstring parameter rename is a **new `lsp_rename` matcher**, not a reconcile
feature. But it does **not** slot into the current matcher contract unchanged —
that contract has a scoping hole for this case.

### The matcher contract must gain location

`argparse.edits(src, old, new)` gets away with no location because
`add_argument` dest strings are structurally distinctive and matched via the
call node. **Docstring param names are not** — `data`, `path`, `name`, `x` recur
across every function in a file. A matcher that scans the whole source for an
entry named `old` would rewrite the `path:` line in *every* function's docstring
when you rename one function's `path`, and would fire when renaming an unrelated
local that merely shares a param name.

The disambiguator is **position**, which the current augment does not pass. So a
prerequisite: extend the augment to hand matchers the triggering edit's location.
When an LSP renames a parameter, the server's edits (declaration + body uses) all
fall inside the enclosing function, so the docstring matcher can scope by *find
the function whose range contains a server edit, touch only its docstring*.
Contract becomes roughly:

```
M.edits(src, old, new, ranges) -> lsp_rename.Edit[]
```

where `ranges` are the server's edited byte ranges in that file. `argparse`
ignores the new arg (its call-structure match already scopes it); the docstring
matcher uses it to pick the one function. This is the "new plumbing" — small, but
real, and a hard prerequisite for the rename milestone.

With that:

- `lua/lsp_rename/docstring.lua` exposing `M.edits(src, old, new, ranges) ->
  lsp_rename.Edit[]`. Find the function enclosing a server edit, locate the entry
  documenting `old` in *its* docstring (via `lua/docstring/`'s
  `locate_doc_entries`), return an edit over just that name token, preserving the
  description.
- Register it: `matchers.python = { require "lsp_rename.argparse", require
  "lsp_rename.docstring" }` in `lsp_rename/init.lua`.

This split **relocates** the hardest ambiguity rather than removing it: with
renames handled at rename-time by `lsp_rename`, the orphan entries reconcile
encounters are (barring a rename done outside LSP rename) genuine *removals*. The
rename-vs-removal ambiguity still exists — it now lives at the lsp_rename
boundary, resolved by the rename *action* supplying the intent that reconcile
alone can't recover.

## Extensibility

Two registries in `lua/docstring/`, keyed independently:

- **language** (`lang/<language>.lua`) → `{ locate_target, signature_params,
  locate_docstring }` (treesitter). Adding a language = new queries.
- **style** (`style/<style>.lua`) → `{ parse, render, locate_doc_entries }` for a
  convention (google / numpy / rest). Adding a style = one reverse parser + one
  renderer (+ the entry locator, which `parse` can expose).

A concrete operation is a `(language, style)` pair. The style default is
configured; auto-detecting an existing docstring's style is later work (attempt
configured style, warn on parse failure).

Naming keeps the two directions distinct: `signature_params` (from code) vs
`locate_doc_entries` (from prose) — never both called "params".

The reverse parsers are the real effort, one per style:

- **Google** — `Args:` section, `name (type): desc`, continuation by indent.
- **NumPy** — `Parameters\n----------`, `name : type` then indented desc.
- **reST** — `:param name:` / `:type name:` lines.

## Testing

Mirror codedocs' golden-file approach (see this repo's `tests/README.md` for the
harness):

- **Round-trip**: `parse(render(m)) == m` per style — the stability guarantee.
  Must include a docstring that *already contains an orphan marker*, so the
  marker is proven to re-parse into `unmatched` rather than corrupt on a second
  pass.
- **Idempotency**: `reconcile(reconcile(x)) == reconcile(x)` — assert the second
  run is a no-op (not that correct docstrings are never touched).
- **Reconcile cases**: add-at-end, add-in-middle, reorder, no-existing-docstring
  (fresh), orphan (preserved + flagged), untyped signature (docstring type
  preserved), `self`/`*args`/`**kwargs`/keyword-only-`*` filtering.
- **Golden fixtures** per `(language, style)`, input→output pairs.

## Integration in this config

- **`lua/docstring/`** — core library: `model`, `style/<style>`, `lang/<language>`,
  `reconcile`, `command`. Reuse `lua/utils/` (cursor/range helpers) rather than
  re-rolling.
- **Command** registered lazily in core (require on first use): `:Docstring`.
- **`lua/lsp_rename/docstring.lua`** — the rename matcher, registered in
  `lsp_rename/init.lua`.

## Milestones

1. **`lua/docstring/` Python+Google**: `DocModel`, `parse`, `render`,
   `locate_doc_entries`, treesitter locators. Round-trip tests. (Foundation for
   both consumers.)
2. **Reconcile Python+Google, end to end**: locate → extract → parse → reconcile
   → render → replace, behind `:Docstring`. Covers add, reorder, fresh-gen,
   idempotency. Full tests.
3. **`lsp_rename` docstring matcher**: rename a param → docstring name token
   tracks it in the same undo block. Reuses the milestone-1 `locate_doc_entries`;
   depends on the augment-contract extension (pass edit ranges).
4. Orphan UX — flag unmatched entries clearly.
5. Second **style** (NumPy or reST) — validates the style boundary for both
   consumers at once.
6. Second **language** — validates the treesitter boundary.
7. Nice-to-haves: tabstops on new blanks; `returns`/`raises` reconcile; class
   docstrings; style auto-detection.

## Open questions

- Orphan rendering: inline marker vs a separate holding section vs virtual-text
  warning + no buffer change? (Marker must round-trip — argues against free-text.)
- Style mismatch on parse failure: fall back to fresh-generate (data loss risk)
  or abort with a warning (safe default)? Leaning abort-and-warn.
- Do we ever touch `summary`/`other` sections, or treat them as strictly
  immutable passthrough? (Default: immutable.)
