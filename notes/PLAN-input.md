# DRAFT PLAN: a rename-specific input float (keep snacks for `vim.ui.input`)

Snacks stays our `vim.ui.input`. We are **not** replacing it globally. The only
consumer whose input is *a term under the cursor* — where char offset, word-start
anchoring and normal-mode start all matter — is **LSP rename**. Rename is the one
case that fits "edit a prefilled term over a buffer", so we route only `grn` to a
small dedicated float and leave every other `vim.ui.input` caller on snacks.

The float is written from scratch — we clearly want something different from what
`snacks.input` is. `snacks/input.lua` is kept only as a **reference to crib the
fiddly, easy-to-botch bits from** (float lifecycle, screen-edge clamping,
winhighlight groups) — not a design dependency, not a fork.

## What we're retiring

The rename work-arounds we ship *inside the snacks input config* (`lua/plugins/pickers.lua`).
Snacks input itself stays enabled and untouched for generic use; only the
rename-specific `win`/`on_win`/`col` block and the `cursor_word_offset` helper go
away — they exist solely to bend snacks into behaving like a rename float, and each
is a work-around for a missing/wrong behaviour in `snacks.input`:

| Want | Work-around we ship | Why it's a hack |
|------|--------------------|-----------------|
| Start in **normal** mode | `on_win` schedules `stopinsert` | `input.lua:256` does an unconditional `startinsert!`; no `mode` option |
| Cursor at the **same char offset** within the prefilled term | `cursor_word_offset` + `place()` re-asserted from a **one-shot `TextChanged`** autocmd | `vim.ui.input` carries no cursor position; and snacks' `expand` handler resets the cursor (below) |
| Float anchored to the **word start**, not the cursor | `col` is a function returning `-2 - cursor_word_offset(current_win)` | `relative="cursor"` anchors to the cursor screen column; word-start anchoring isn't an option |

The re-assert exists purely to undo snacks' own damage: on the first `TextChanged`
(inserting the default text), the `expand` handler (`input.lua:268-284`) calls
`win:update()` → `nvim_win_set_config` (`win.lua:803`), which **snaps the buffer
cursor to column 0**. We register a `TextChanged` autocmd *before* snacks' and
schedule the cursor set to run *after* its reset. That ordering dependency is the
most brittle part.

We also lean on `winnr("#")` to recover the parent window inside `on_win` (the float
already has focus by then) — snacks captures `parent_win` at `input.lua:104` but
never exposes it.

## Root causes (why we don't build rename on top of snacks input)

All of these live in `snacks` and can't be configured away, so a dedicated float is
cleaner than fighting them:

1. **Unconditional `startinsert!`** — no way to open in normal mode.
2. **`prompt` buftype** — forces the insert-mode cursor to the prompt end and is why
   normal-mode cursor placement is fragile. The `prompt_setcallback`/`prompt_setinterrupt`
   machinery exists only because the buffer is a prompt buffer.
3. **`expand` resets the cursor** — `win:update()` fully reconfigures the float
   (recomputing `col`) on every text change, losing cursor/view state.
4. **No cursor-offset in the API** — `vim.ui.input` opts have `default` but no
   "place cursor at column N", so rename can't say "cursor was mid-word".
5. **`parent_win` not exposed** — hooks can't get the origin window.

## Scope: rename only, wired at the `grn` keymap

`grn` is the built-in mapping (`vim/_core/defaults.lua:204`): `vim.lsp.buf.rename()`
with no arg, which prompts via `vim.ui.input` (→ snacks). We override `grn` to open
our float instead and hand the confirmed name straight to `vim.lsp.buf.rename(new_name)`,
which **skips the `vim.ui.input` prompt entirely** when `new_name` is supplied
(`buf.lua:787` / `:841`). So snacks input is never involved in rename, and no other
consumer changes.

`vim.lsp.buf.rename` captures its origin window early (`win = nvim_get_current_win()`,
`buf.lua:736`) and derives the position params from it, so our float stealing focus
is fine — the rename still targets the original symbol.

- **Every other `vim.ui.input`** (`<leader>:!`, agentic.nvim, nvim-tree, mini.files,
  neogit, nvim-dap, diffview, metals, flutter-tools, fzf-lua, gitsigns, …) stays on
  snacks, unchanged.
- **Rename** uses the new float: normal mode, cursor at the char offset within the
  term, float anchored to the word start.

Because rename is invoked directly (not via a global `vim.ui.input` shim), there is
**no need to detect "is this a rename?"** — no `default == cword` self-gating. The
call site *is* the gate.

**Trade-off:** passing `new_name` to `rename()` takes the early-return branch
(`buf.lua:787`), skipping the `LspReferenceTarget` highlight that the interactive
path draws over the symbol (`buf.lua:806`). We don't draw that highlight today
(snacks path), so this is no regression; add it in the float later if wanted. The
`prepareRename` validity check still runs — the `new_name` return happens *inside*
the `prepareRename` callback — so invalid rename targets are still rejected.

## The one decision that dissolves most of this

**Use a one-line modifiable scratch buffer, not a `prompt` buffer.** With a normal
buffer we own the mode and cursor outright: no `startinsert!`, no prompt-end
snapping, no prompt callbacks, no `expand`-driven reconfigure. Confirm/cancel become
plain buffer-local keymaps instead of prompt callbacks.

## Design — `modules/wordinput.nvim`

A focused float for editing the term under the cursor. The sole caller today is the
`grn` remap, but the API is generic enough to reuse for any future cword-anchored
input.

`require("wordinput").open(opts, on_confirm)` — signature mirrors `vim.ui.input`
(`opts.prompt`, `opts.default`), so the rename call site reads naturally. It assumes
`opts.default` is the term at the parent cursor (true for rename) and derives the
offset/anchor from the live cursor position.

**Written directly, not adapted from snacks.** Most of it is trivial: a scratch
buffer, one `nvim_open_win`, `title = opts.prompt` (no icon), `border`/`winhighlight`
set on the win config, two buffer-local maps. Nothing to lift there.

**Consult `snacks/input.lua` as a reference only for the fiddly parts** — where
getting the details right is the hard bit and cribbing beats reinventing:
- **Window lifecycle / teardown** — clean close, autocmd + keymap cleanup, guarding
  double-close, and treating focus-loss (`WinLeave`/`WinClosed`) as an implicit
  cancel so a stray click doesn't strand the float.
- **Screen-edge clamping** of the `relative="cursor"` placement, so the float doesn't
  spill off-screen near the bottom/right edge.
- The exact **winhighlight** group mapping.

**`restore-parent-window-and-mode-on-close` — verify the scenario first, else omit.**
This is not a code-quality question about snacks; it's whether our *normal-mode plain
float* ever leaves the parent in a wrong state. `nvim_win_close` should return focus
and cursor to the parent cleanly, so there is likely nothing to restore. Confirm
that; only add restore logic if a concrete wrong-state case shows up.

**Drop entirely:** the prompt icon, prompt buffer + `prompt_set*`, `startinsert!`,
`expand`/`win:update`, history, completion (`complete`/`completefunc`), backdrop.

`M.open(opts, on_confirm)`:

1. **Capture origin up front** (before any window work — parent is unambiguously
   current here, so no `winnr("#")`): `parent_win`, and the `<cword>` + its start
   column at the parent cursor (lift the verified `\k*$` / `iskeyword` match out of
   `pickers.lua` into the module).
2. **Compute offset/anchor**: `offset` = cursor's byte offset within the word;
   `anchor` = word-start column.
3. **Buffer**: scratch, `modifiable`, one line = `opts.default or ""`.
4. **Float**: `relative="cursor"`, `col = base_col - (cursor_col - anchor)`, `row`
   above the cursor, fixed `width`.
5. **Cursor + mode**: `nvim_win_set_cursor(win, {1, offset})`; stay in normal mode.
   Plain buffer ⇒ one synchronous call, no `vim.schedule`, no re-assert.
6. **Keymaps** (buffer-local): `<CR>` in **both normal and insert** mode (an
   insert-mode `<CR>` would otherwise split the one-line buffer) → read line, close,
   `on_confirm(text)`; `<Esc>`/`<C-c>` → close, `on_confirm(nil)`. **Close before
   `on_confirm`** so `nvim_get_current_win()` inside the ensuing `rename` resolves to
   the parent.
7. **Cleanup**: close the window on confirm/cancel, guarded by `pcall`/validity checks.

### Explicitly NOT in scope

Completion, history, and backdrop (a dimming overlay behind the float — snacks
already disables it for input). Add any of these later only if a real need
appears; each would be a clean addition, not a work-around.

### Follow-up: size the float to the term (grow/shrink)

Currently the float is a fixed `width=30`; a term longer than the box scrolls the
view instead of widening. Size the width to the edited text and adapt it as the
term grows/shrinks — the `expand` feature snacks had, but clean here because the
plain buffer removes the machinery that made snacks' version reset the cursor.

Verified feasibility (headless, nvim 0.12):
- A **width-only** `nvim_win_set_config(win, { width = W })` on the *focused*
  float leaves its position unchanged (does **not** re-anchor to the float's own
  cursor — stays at the word-start `col` set at open) and preserves the buffer
  cursor. So a `TextChanged`/`TextChangedI` autocmd recomputing
  `W = vim.api.nvim_strwidth(line) + 1` (the `+1` keeps a cell for the cursor when
  appending) is sufficient. No `win:update`-style full reconfigure ⇒ none of
  snacks' cursor-reset applies.
- **Right-edge alignment is *why* to size to the term, not a separate concern.**
  `:help api-floatwin` (`api.txt:3844`): the builtin UI truncates a float's
  position so it's fully within the screen grid — i.e. it shifts the float left to
  fit. For the current fixed `width=30` float that shift breaks word-start
  alignment whenever the cursor is within 30 cols of the right edge (the reported
  bug). A float sized to the term is narrow, so it only gets shifted when the term
  *itself* would overflow the screen edge — unavoidable, and rare. So term-sizing
  is the fix for the misalignment; there is nothing extra to clamp manually.

Set the initial `width` from `opts.default` by the same formula so open and resize
agree.

### Why this removes every work-around

- Normal-mode start → native (no `startinsert!`).
- Cursor at offset → set once, synchronously; nothing resets it.
- Word-start float anchor → `col` from the captured offset at open.
- Parent window → captured directly; `winnr("#")` gone.
- The one-shot `TextChanged` re-assert and its registration-order dependency → gone.

## Migration steps

1. Create `modules/wordinput.nvim/lua/wordinput/init.lua` per the design; lift
   `cursor_word_offset` out of `pickers.lua` into the module.
2. Add the module to `rtp` in `init.lua` (as other `modules/*` dev plugins). No
   `setup()` needed — it exposes `open()` and nothing global.
3. In `pickers.lua`: **leave `input` enabled** (snacks stays our `vim.ui.input`).
   Delete only the rename-specific `win`/`on_win`/`col` block and the
   `cursor_word_offset` helper. Keep any generic styling (icon, width) as taste.
4. Remap `grn` (in `plugin/keymaps.lua`, near the other `gr*` maps) to:
   ```lua
   local wi = require("wordinput")
   local cword = vim.fn.expand("<cword>")
   wi.open({ prompt = "New Name: ", default = cword }, function(name)
       if name and #name ~= 0 then vim.lsp.buf.rename(name) end
   end)
   ```
5. Verify:
   - **Rename** (the cases we debugged): normal mode; cursor lands at the word offset
     and **stays** (the old bug); float anchored to word start for cursor at char 1,
     2, 3…; `<CR>` renames via LSP; cancel via `<Esc>` leaves the buffer untouched.
   - **Generic unchanged**: `<leader>:!` and a plugin prompt (nvim-tree create/rename)
     still use snacks input exactly as before.
6. Once stable, update `CLAUDE.md` (drop the snacks rename-input work-around notes)
   and remove this plan.

## Testing notes (learned the hard way)

- Since `grn` invokes the float directly (not via `vim.ui.input`), it's headless-
  drivable without the snacks UI-event gating that blocked the old path. But the
  **cursor-reset bug only reproduced with a live UI** (the offending
  `TextChanged`→reconfigure didn't fire headless), so do a real-UI check (kitty pty
  or manual) for final confirmation.
- Byte-offset semantics: `nvim_win_set_cursor` wants byte columns and `\k*$` is
  byte-length — consistent for ASCII identifiers (the common case).
