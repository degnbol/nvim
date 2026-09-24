# Plan: a documentation tier in the K chain
Status: spec

## Relevant files

- `plugin/keymaps.lua:428-503` — `peek_file`, `lsp_hover_capable` and the `K`
  chain the new tier joins.
- `lua/utils/init.lua:399` — `cword_start_col`, `match_covering`.
- `ftplugin/markdown.lua` — its comment records that core gives hover floats
  `filetype=markdown`.
- `/opt/homebrew/share/nvim/runtime/lua/vim/lsp/util.lua:1686-1810` —
  `open_floating_preview`: anchor buffer, `lsp_floating_bufnr`, close autocmds.
- `~/.claude/skills/neovim/references/ui.md` (float lifecycle),
  `references/keymaps.md:38` (the only description of the `K` chain, already
  stale).
- `.claude/skills/lsp/hover_probe.py` — for checking a server's raw reply while
  building this.
- `tests/README.md`.
- `notes/PLAN-stub-docstrings.md` — companion plan.

## Background

With `PLAN-stub-docstrings.md` landed, hover on `counts.sum` ends at a pointer:

```
Refer to `numpy.sum` for full documentation.
```

Resolving that inside the stub would take `numpy/__init__.pyi` to ~5× pristine,
because the text multiplies by overload count. Resolving it on demand costs
nothing.

`KK` already reaches the float — `vim.lsp.buf.hover` passes
`focus_id = "textDocument/hover"`. `K` on `numpy.sum` inside it does not work:
the float has no LSP client, so `lsp_hover_capable(0)` is false and the chain
falls through to `keywordprg`. Measured inside a hover float under real startup:

```
filetype             = markdown          (not empty: core stylizes when syntax is on)
lsp clients          = 0
buftype              = nofile
w.lsp_floating_bufnr = <source bufnr>
keywordprg           = :Man
<cword> on numpy.sum = sum
```

### The lookup reuses the running client

A dotted name has no position in any document, but it can be given one: a
throwaway buffer holding `import numpy` and `numpy.sum`, attached to the
basedpyright client that is already serving the source buffer, answers a normal
`textDocument/hover`. Verified against the live client — nothing new is started,
no interpreter has to be resolved, and no subprocess runs. `bufadd` + `bufload`
never touch disk, so the scratch path needs no cleanup.

A synthetic reference is a *bare* reference, so the reply carries every overload:
12048 chars over 359 lines for `numpy.sum`. Adding a call form does not help —
with no arguments to match, no overload wins. basedpyright separates signatures
from prose with a `---` line, and keeping only what follows it gives 4233 chars
over 116 lines, starting `Sum of array elements over a given axis.` That is the
right half regardless: the signature was already on screen in the hover the
pointer came from.

### Other facts the design rests on

- **The originating buffer is already recorded.** `open_floating_preview` sets
  `vim.w[win].lsp_floating_bufnr` unconditionally for every float it creates
  (`util.lua:1776`), unlike the `focus_id` variable, which exists only when one
  was passed. So the float's language and its client are both recoverable with no
  new state, and it works for peek and signature-help floats too.
- **The float's own filetype is `markdown`**, a plausible-looking wrong answer.
  The origin buffer must outrank it.
- **`<cword>` is not enough.** `.` is not in `'iskeyword'` (measured
  `@,48-57,_,192-255`), so `numpy.sum` yields `numpy` or `sum` depending on the
  cursor. numpy also renders `See Also` entries unquoted —
  `numpy.sum : equivalent function` — so the reader must not rely on a code span.

## Changes

### `lua/utils/` — three general helpers

None carries documentation-lookup knowledge, so none belongs in the new module:

```lua
-- lua/utils/init.lua, beside cword_start_col
--- Dotted identifier around the cursor, `numpy.sum` where <cword> gives `sum`.
--- The first segment must start `[%a_]`, so `np.float64(1.5)` does not read
--- `1.5`; leading and trailing dots are dropped.
--- @param winid integer 0 = current
--- @return string|nil
function M.dotted_cword(winid) end

-- lua/utils/lsp.lua
--- Buffer an LSP floating preview was opened from, nil outside one.
--- @param winid integer 0 = current
--- @return integer|nil bufnr validated with nvim_buf_is_valid
function M.float_origin_buf(winid) end

-- lua/utils/treesitter.lua
--- @param winid integer 0 = current
--- @return string|nil lang
function M.language_at_cursor(winid) end
```

### New `lua/refdoc.lua`

The scratch buffer, the request and the dispatch. (Not `doc_lookup`:
`lua/docstring/` already exists and edits docstrings in-buffer.)

```lua
--- Documentation for a dotted name, from a client already serving the language.
---
--- Never starts a client: a language with none running has no answer here.
--- @param name string dotted symbol
--- @param client vim.lsp.Client
--- @param on_result fun(lines: string[]|nil, err: string|nil)
function M.request(name, client, on_result) end

--- Look up the symbol under the cursor against the language it belongs to.
--- The buffer that answers is a float's origin buffer where there is one, else
--- the window's own. Returns synchronously whether it dispatched, so the caller
--- can fall through to its next tier without waiting on the callback.
--- @param winid integer 0 = current
--- @param on_result fun(lines: string[]|nil, err: string|nil)
--- @return boolean dispatched
function M.lookup(winid, on_result) end
```

- The scratch buffer takes a path inside the client's root so the URI is one the
  server accepts, holds `import <root segment>` and the name, and is attached
  with `vim.lsp.buf_attach_client`. Reuse one buffer per client rather than
  creating one per keypress.
- Keep only the text after the `---` separator; where there is none, keep all of
  it.
- Cap the result and pass `max_width`/`max_height`, as `peek_file` does with
  `PEEK_LINES`.
- Never silently swallow: `on_result(nil, err)` ⇒ `vim.notify(err, WARN)`, with
  distinct text for "nothing resolved that name" and "the request failed".
- Python only at first; a second language is one more branch. Normalise compound
  filetypes — `python.blender` is configured in `after/lsp/basedpyright.lua:83`.

### `plugin/keymaps.lua`

- `DOC_FOCUS = "refdoc"` beside `PEEK_FOCUS`; the `wincmd p` branch at the top of
  the `K` callback tests both. The reason is not the round trip —
  `open_floating_preview` bounces focus itself (`util.lua:1696-1702`) — it is to
  return before issuing a request.
- The tier sits **after** `lsp_hover_capable`, except where the treesitter
  language at the cursor differs from the buffer's filetype — an injected region,
  such as a ```python fence in a markdown note. Without that exception the tier
  never runs there, because marksman answers `textDocument/hover` for markdown
  and wins the branch above. The exception costs nothing elsewhere, since the two
  agree in an ordinary buffer.
- The callback opens `open_floating_preview(lines, "", { focus_id = DOC_FOCUS })`
  — empty filetype, since markdown normalisation would eat the docstring's
  underlines and indentation.
- The reply arrives outside the keypress, so `vim.schedule` it and re-check
  before opening: the current window is still `winid`, the cursor has not moved,
  and `dotted_cword(winid)` is still the requested name. Without this,
  `open_floating_preview` takes its anchor from whatever buffer is current at
  call time (`util.lua:1686`) and closes that buffer's `lsp_floating_preview`
  (`:1716-1719`), so a late result would close an unrelated float. Drop an older
  in-flight lookup when a new one starts.

Two floats coexisting is safe, verified: opening from inside the hover float
reads `vim.b[float_buf].lsp_floating_preview` (nil), so the hover float survives,
and the doc float is closed by the hover float's own autocmds.

### Tests

Per `tests/README.md`, without a live LSP: `dotted_cword` on `` `numpy.sum` ``,
on the bare `numpy.sum : equivalent function` form, at every character of the
name, on a trailing dot, and on `np.float64(1.5)`; `float_origin_buf` preferring
the recorded bufnr over the float's own `markdown` filetype, and returning nil for
a wiped origin; the `---` split keeping the prose half and tolerating its absence.

### Documentation

`~/.claude/skills/neovim/references/keymaps.md:38` is the only description of the
`K` chain and already cites a stale line number. Correct it with the change,
naming the tiers rather than counting them. It lives in the global skills tree,
not this repo, so it is a separate commit.

## Expected outcome

`KK` into the hover float on `counts.sum`, `K` on `numpy.sum`, and the
`numpy.sum` documentation opens in a second float, 116 lines of prose with no
signature pile; `K` there returns focus. Zero stub-size cost, and nothing
launched that was not already running.

The tier also fires on a dotted name in a Python comment or docstring wherever no
client answers hover, and inside a Python fence in a markdown note — using the
client already serving some Python buffer in the session.

## Non-goals

- **Starting a language server to answer a lookup.** A language with no client
  running has no answer here. That is what keeps markdown injections cheap: at
  most one extra buffer on a client that already exists.
- **Languages beyond Python at first.** Lua already has `ftplugin/lua.lua:25`.
- **Changing when LSP hover wins**, beyond the injected-region exception above.
- **Chasing a name from inside the doc float** — `K` there returns focus, as in
  the peek float.
- **Rendering the result as markdown.**

## Open questions

None.
