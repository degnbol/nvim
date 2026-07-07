# Plan: live line reflow for wordinput.nvim

## Goal

`wordinput.nvim` opens a transparent float over the word under the cursor as an
input surface (used for LSP rename, etc.). The float overlays the word in place
and grows/shrinks as you type. Today the **trailing text on the line does not
move**: when the typed term grows past the original word, the float's opaque
paint covers the code after the word (`= computeValue()` etc.); when it shrinks,
the original word's tail peeks past the float's right edge.

Make the trailing buffer text on the word's line reflow live so it always begins
immediately after the float, whether the term is longer or shorter than the
original word.

File to edit: `lua/wordinput/init.lua` (single file, ~184 lines). Read it first.

## Hard constraints (do not violate)

- **Never mutate the real buffer.** The buffer must stay pristine until the
  actual LSP rename fires with `(position, newWord)` on the original text.
  Reflow must be display-only (extmarks in the parent buffer).
  - Rationale: LSP document sync attaches via `nvim_buf_attach(bufnr, false,
    {on_lines=…})` (see `runtime/lua/vim/lsp.lua`), which is *below* the autocmd
    layer. `noautocmd`/`eventignore` and low-level `nvim_buf_set_text` do **not**
    suppress it — every edit notifies the server, thrashing diagnostics/semantic
    tokens against half-typed identifiers. So the "edit-and-revert" approach was
    rejected. Do not revive it.
- **Never change window options** (`conceallevel`, `concealcursor`, etc.) on the
  parent window. Forcing `conceallevel` would perturb conceal of unrelated
  content across that window's whole viewport. The reflow must respect whatever
  the user's settings already are, and degrade gracefully.
- **No settings inspection / no gate.** The conceal extmark below is self-gating
  (a no-op at `conceallevel = 0`), so there is no need to read `conceallevel` at
  all. Set the extmark unconditionally.

## Why display primitives, and the width levers

To reflow without editing the buffer you must change the *rendered* width of the
region the word occupies. Neovim has exactly two display levers for that:

- **inline `virt_text`** (`virt_text_pos = "inline"`) — *adds* rendered width.
- **`conceal = ""` extmark** — *removes* rendered width, but only when the
  window's `conceallevel >= 1`; at `0` it renders normally (this is what makes it
  self-gating).

There is no third lever and no way to remove width independent of
`conceallevel`. This drives the asymmetric design below.

## The design

Let `orig_width = vim.api.nvim_strwidth(default)` (display width of the original
word), captured once at open. On every keystroke compute
`new_width = vim.api.nvim_strwidth(current_line)` and `delta = new_width -
orig_width`, then maintain **one** extmark (a single stable `id`) in the
**parent** buffer, reconfigured by the sign of `delta`:

- `delta > 0` (term longer): inline `virt_text` of `delta` spaces at the word's
  end column (`word_start + orig_word_byte_len`). Trailing text is pushed right
  by `delta` so it begins where the grown float ends. **Always works** — no
  dependency on any setting.
- `delta == 0`: clear the extmark (delete it).
- `delta < 0` (term shorter): conceal the last `-delta` cells of the original
  word — an extmark spanning `[word_end - (-delta), word_end)` with
  `conceal = ""`. This collapses those cells so trailing text slides left to meet
  the float. Behaviour by the user's existing `conceallevel`:
  - `>= 2`: cells fully hidden → exact.
  - `== 1`: cells collapse to a single replacement cell → trailing lands one cell
    too far right (a faint `Conceal`-highlighted space). Accepted as cosmetic.
  - `== 0`: no-op → original tail peeks past the float, i.e. today's behaviour.
    Graceful, no surprise.

Note the shrink path works in *display cells* but conceal takes a *byte* range;
map the cell count to the corresponding byte range of the original word's tail
(the word may contain multibyte chars). For the common ASCII-identifier case
cells == bytes, but do it correctly.

## Where to hook

The resize logic already lives in two autocmds on the float buffer — reuse them,
don't add new ones:

- `InsertCharPre` (init.lua ~113): fires before `v:char` lands; currently widens
  the float ahead of the char. The reflow's grow case wants the same
  look-ahead so trailing text moves in the same frame the float widens — compute
  `delta` against `line .. vim.v.char` here for the grow pad.
- `TextChanged` / `TextChangedI` (init.lua ~123): authoritative resize+repaint
  for deletes and normal-mode edits. Recompute `delta` from the actual buffer
  line and reconfigure the extmark here (this is where shrink is handled;
  `InsertCharPre` never fires on deletion).

Factor a single `reflow(width)` helper (mirroring the existing `set_width` /
`paint` helpers) that both handlers call, so grow/shrink/clear logic lives in one
place.

## Anchoring details

- The parent buffer, parent window, `word_start` (byte col of word start), and
  the original word text are all available at open (init.lua ~48–54). Capture
  `orig_width` and the word's end byte column there too.
- Reflow extmarks go in the **parent** buffer, in a namespace distinct from the
  float's existing `"wordinput"` paint namespace (or reuse it with a separate
  `id` — but a separate namespace is cleaner for teardown).
- **Teardown:** delete the reflow extmark(s) when the float closes. The existing
  `finish()` / `WinClosed` path (init.lua ~136–162) is where cleanup belongs. The
  parent buffer outlives the float, so a leaked extmark would visibly corrupt the
  line — make teardown unconditional (also on the `WinLeave` cancel path,
  init.lua ~164).

## Verify empirically (scratch buffer, do not guess)

1. **Grow:** rename a short word to a longer one on a line with trailing code;
   confirm trailing code slides right and stays visible as you type, and returns
   to place on cancel.
2. **Shrink at `conceallevel = 2`:** confirm trailing code slides left exactly to
   the float edge.
3. **Shrink at `conceallevel = 0`:** confirm no crash and the tail simply peeks
   (degradation), extmark cleaned up on close.
4. **Shrink at `conceallevel = 1`:** confirm the 1-cell residual and check
   whether it lands under the float's opaque term paint (invisible) or aligns
   with the float's transparent trailing caret cell (shows through as a faint
   space). Report which; if it peeks and is distracting, adjust the concealed
   cell count by one.
5. **Non-current-window conceal:** the word is on the *parent* window's cursor
   line while the float holds focus. Confirm conceal actually renders there —
   `'concealcursor'`'s "reveal on cursor line" is documented for the *active*
   window, and the parent is not active, so it should conceal normally, but
   verify. If it does not conceal, the fallback is the same graceful degradation
   (tail peeks), so correctness never breaks — only the shrink bonus is lost.
6. **Multibyte:** rename involving a multibyte original word to confirm the
   cell→byte mapping for the conceal range is correct.
7. Cancel (`<Esc>`/`<C-c>`), confirm (`<CR>`), and focus-loss (`WinLeave`) each
   leave the parent line pristine with no stray extmark.

## Out of scope

- Reflowing anything but the single line the word is on.
- Syntax-highlighting the pushed/collapsed text differently.
- Any change to the LSP rename call itself — it still runs on the untouched
  buffer via the existing `on_confirm(value, col)` contract.
