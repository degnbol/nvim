# Plan: `K` peeks file contents under cursor

## Goal

Overload normal-mode `K` so that, when the cursor is on a path that resolves to
a readable file, it opens a floating preview of that file — the same way `K`
already peeks info on a symbol via LSP hover. Primary use case: a path in a
markdown doc, a comment, a log, a shell/config file, where LSP hover has
nothing useful to offer.

## Background / current state

- `K` is the **neovim built-in default**, not a custom mapping. With an LSP
  attached, neovim installs a buffer-local `K` → `vim.lsp.buf.hover()` — but
  **only if no `K` mapping already exists** (`runtime/lua/vim/lsp.lua:868-873`,
  `maparg('K','n') == ''` check). Without an LSP, `K` runs `keywordprg`.
- The custom hover-ish keymap `gh` → `vim.lsp.buf.signature_help`
  (`lua/autocmds/lsp.lua:57`) is unrelated.
- **blink.cmp is irrelevant** — it is completion only and never touches `K`/hover.
- There is already a stub for this intent: `plugin/keymaps.lua:390-412` has a
  commented-out `LspAttach` `K` map plus a TODO ("the LSP nmap K should also
  show the float aligned to the function name start and not cursor"). That stub
  was buffer-local — wrong for this feature (see Decision 7).
- There is already path-resolution logic to reuse: the custom `gf` in
  `plugin/paths.lua:79-126` extracts a path token under the cursor and expands
  `$(git root)`, `$VAR`, and buffer-local `NAME=value` assignments, falling back
  to `normal! gf`. This is the resolver this feature needs.

## Key constraint discovered

`vim.lsp.buf.hover()` has **no fallback hook**. It is async
(`buf_request_all` + callback); on empty results it just calls
`vim.notify('No information available')` and returns — no return value, no
`on_empty` config, no continuation callback (`runtime/lua/vim/lsp/buf.lua:75`+).
So a true "try hover, and *if empty* peek the file" is not cleanly available.
The decision must be made **synchronously, up front**, before calling hover.

## Decisions (from grilling)

1. **Fallback, not a new key.** Overload `K`; do not introduce a separate
   keymap. `K` already means "tell me about the thing under the cursor."

2. **File-first discriminator.** Decide synchronously: if the token under the
   cursor resolves to a readable file → peek it; otherwise → hover (then
   keywordprg). LSP symbols essentially never resolve to readable file paths, so
   the overlap risk is negligible. No async empty-detection trickery.

3. **Resolution chain.** Turn the cursor token into a path, first readable hit
   wins:
   1. `$(git root)` / `$VAR` / buffer-`NAME=value` expansion (lift from the
      existing `gf` in `plugin/paths.lua`).
   2. `expand('<cfile>')` with `~`/env expanded; if absolute & readable, use it.
   3. Join with the current buffer's directory.
   4. `vim.fn.findfile(cfile)` to honour `'path'` / `'suffixesadd'`.
   Buffer-relative-first matches the doc/comment use case.

4. **Performance — non-issue.** Worst case is `findfile`'s `'path'` walk. The
   user's `'path'` is a flat list of ~40 dirs (mostly leaked `$PATH` bins) with
   **no recursive `**`** — only single-level `./src` / `src/*`. So ~40 `stat()`
   calls worst case, sub-millisecond on local SSD, and `K` is a deliberate
   keypress, not a per-keystroke hot path.

5. **Cheap token-shape gate.** Skip the filesystem walk entirely for bare
   identifiers (the common LSP-hover case):
   ```lua
   if not cfile:match('[/~]') and not cfile:match('%.%w+$') then return hover() end
   ```
   `myFunction` costs one string match and zero stats. `foo.bar` (method access)
   matches `%.%w+$` and attempts resolution, but won't resolve to a readable
   file, so it falls through to hover after one cheap `filereadable` miss.
   Acceptable; flip the order later if it ever feels slow.

6. **Render with `vim.lsp.util.open_floating_preview()`.** Already in the
   runtime, already the machinery hover uses — same border, max width/height,
   `q`/`<Esc>` close, focus-on-repeat. It has **no `ft` opt** — the signature is
   `(contents, syntax, opts)`, where `syntax` sets `vim.bo.syntax`. Pass the
   matched filetype (`vim.filetype.match({ filename = path })`) as that `syntax`
   arg for regex highlighting. Treesitter would need an extra
   `vim.treesitter.start` — deferred (YAGNI).

7. **Content cap + anchors.**
   - Read at most the first 500 lines: `vim.fn.readfile(path, '', 500)` —
     bounded cost even for a huge file, more than fills the float.
   - Strip a trailing line/anchor suffix (`:<lnum>` or `#anchor`) **during
     resolution** so `foo.lua` followed by `:42`, and `README.md#install`, still
     pass `filereadable`. **Defer** scrolling the peek to that line (YAGNI — add
     cursor positioning when actually missed).

8. **Single global `K`.** Define unconditionally in `plugin/keymaps.lua`,
   replacing the dead stub at 402-412. Must be **global, not buffer-local**:
   - A no-LSP buffer (markdown/log) never fires `LspAttach`, so a
     buffer-local-on-attach `K` would not exist there — but that is the primary
     use case.
   - Defining it at startup means the `LspAttach` default sees a `K` mapping and
     skips installing its own, so our mapping stays authoritative everywhere and
     *we* call hover as the fallback.

9. **Built-in `K` via noremap feedkeys.** The keywordprg fallback (case 3) is
   preserved with no recursion:
   ```lua
   vim.api.nvim_feedkeys(vim.v.count1 .. 'K', 'nx', false)  -- n = noremap, x = execute now
   ```
   The `n` flag bypasses all mappings (global + buffer-local), invoking the
   built-in `K`.

10. **Fold in the existing alignment TODO.** The hover branch incorporates the
    `offset_x = cword_start - c` alignment the 390-stub already wanted (align the
    hover float to the word start, not the cursor column).

## Implementation

### 1. New `utils` function — `resolve_path_under_cursor(bufnr) → string|nil`

Lives in `lua/utils/` (genuinely two callers — `gf` and `K` — so it clears the
"cross-file reuse" bar in CLAUDE.md; pick an existing module or a small new one).
Pure: returns an absolute readable path or `nil`. Steps:

1. Token-shape gate (Decision 5) — return `nil` fast for bare identifiers.
2. `$(git root)` / `$VAR` / buffer-`NAME=value` expansion (lifted from
   `plugin/paths.lua` `gf` body).
3. Plain `<cfile>` chain (Decision 3): `~`/env expand → buffer-dir-relative →
   `findfile`.
4. Strip trailing line/anchor suffix before the `filereadable` test; discard
   the anchor for v1.
5. `nil` if nothing readable.

### 2. `gf` (`plugin/paths.lua`) — rewire to the helper

```lua
local p = require("utils.…").resolve_path_under_cursor(bufnr)
if p then vim.cmd.edit(vim.fn.fnameescape(p)) else vim.cmd("normal! gf") end
```

`normal! gf` stays as the final net (covers `includeexpr` / edge cases the
helper does not). Strictly a superset of today's behaviour — no regression.

### 3. Global `K` (`plugin/keymaps.lua`) — replace the stub at 402-412

```lua
vim.keymap.set('n', 'K', function()
  local p = require("utils.…").resolve_path_under_cursor(0)
  if p then peek_file(p) return end                       -- new code, see below
  if hover_capable(0) then
    local cword_start = require("utils.init").cword_cols() -- folds in the 390 TODO
    local c = vim.api.nvim_win_get_cursor(0)[2]
    return vim.lsp.buf.hover({ offset_x = cword_start - c })
  end
  vim.api.nvim_feedkeys(vim.v.count1 .. 'K', 'nx', false)  -- keywordprg, noremap
end, { desc = 'Peek file / hover / keywordprg' })
```

`peek_file(p)` — the only genuinely new code, inline in keymaps.lua (single
caller): strip anchor, `readfile(p, '', 500)`,
`vim.lsp.util.open_floating_preview(lines, '', { ft = vim.filetype.match({ filename = p }) })`.

`hover_capable(0)` — true iff some attached client for buf 0 supports
`textDocument/hover` (`client:supports_method('textDocument/hover')`).

## Deferred (YAGNI)

- Scroll-to-line for the `:<lnum>` / `#anchor` suffix — strip now, jump later.
- Flipping resolve-vs-hover order if the `filereadable` miss ever feels slow.
- Promoting `peek_file` to `utils` — only if a second caller appears.

## Resolved implementation details

- `resolve_path_under_cursor` lives in a new `lua/utils/paths.lua`, alongside
  `git_root` and `buffer_var` (both moved out of `plugin/paths.lua`, which now
  imports them — they have callers in both startup code and the resolver).
- `K` checks `fs_stat(path).type == "file"` before peeking: the resolver also
  returns directories (via the `$(git root)` / `$VAR` case, which `gf` wants),
  so dirs fall through to hover.
