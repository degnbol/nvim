# kitty syntax: vim-kitty vs neovim's builtin

`after/syntax/kitty.vim` is one line, ending `nextgroup=kittyPath`, and no
`kittyPath` syntax item exists — `syntax list kittyPath` is empty,
`nvim_get_hl(0, {name = "kittyPath"})` is `{}`, and `synstack()` over a
`startup_session` path reports `kittySt`. That `nextgroup` has been inert since
the vim.pack migration: the group came from an uncommitted in-place edit to
vim-kitty's own `syntax/kitty.vim` in the lazy.nvim checkout, now deleted. The
lost diff:

```diff
-syn match kittyInclude '^\s*\(env\|glob\)\?include' display
+syn match kittyPath '.*$' contained
+syn match kittyInclude '^\s*\(glob\|gen\)\?include' display nextgroup=kittyPath
+syn match kittyInclude '^\s*envinclude' display

 hi def link kittyNumber Number
+hi def link kittyPath String
```

**Do not port it.** Its motivation is already met upstream, just not reachable.

## Why the patch existed

vim-kitty routes an option's value through `kittySt`, which has no `hi def link`
— so values render as plain text. Only numbers and colours inside them pick up a
group. `kittyPath` → `String` was to make paths visible. The second half split
`envinclude` out of `\(env\|glob\)\?include` so that `geninclude` stopped being
flagged `kittyInvalidKeyword`; it is absent from vim-kitty's generated keyword
list.

## Neovim ships its own kitty syntax

`$VIMRUNTIME/syntax/kitty.vim` (2025-09) highlights *every* option's value
generically via `kittyOption` → `kittyOptionName` (Keyword) + `kittyOptionValue`
(String), so no per-keyword path group is needed. Measured with `--clean`:

| line | name | value |
|---|---|---|
| `startup_session /path/to/x.conf` | `kittyOptionName` | `kittyString` |
| `globinclude glob/*.conf` | `kittyOptionName` | `kittyString` |
| `geninclude gen.py` | `kittyOptionName` | `kittyString` |
| `envinclude KITTY_*` | `kittyOptionName` | `kittyString` |

Both halves of the patch fall out for free, and `nvim --clean` already detects
`kitty.conf` as `ft=kitty`, so vim-kitty is not needed for ftdetect either.

**It is currently shadowed.** Runtime order is vim-kitty → `$VIMRUNTIME` →
`after/`, and vim-kitty sets `b:current_syntax`, so the builtin's
`exists("b:current_syntax")` guard makes it `finish`. The better file never runs.

## The actual decision

Keeping vim-kitty buys two things the builtin lacks:

- `kittyInvalidKeyword` / `kittyInvalidAction` → `Error`, i.e. unknown option and
  action names are flagged. The builtin validates nothing — `not_a_real_option 5`
  highlights as a normal `kittyOptionName`.
- `syntax/kitty-session.vim`, which has no builtin counterpart.

Against that: its keyword list is hand-generated and goes stale (`geninclude` is
already missing), values are unhighlighted, and the builtin's `map` /
`mouse_map` / `kitty_mod` / line-continuation / colour / alpha handling is
considerably richer.

Leaning towards dropping vim-kitty and taking the builtin, accepting the loss of
option-name validation — `modules/kitty-conf.nvim` already carries a current
`kitty_options.json` and is the better place to rebuild validation if it is
missed.

Either way `after/syntax/kitty.vim` can go: `startup_session` is already in
vim-kitty's keyword list, so the line is redundant under vim-kitty and
unnecessary under the builtin.

One builtin wart to expect: in `include ./other.conf` the leading `.` matches
`kittyNumber` (`/[+\-*\/]\{0,1}[0-9.]\+/`), so it colours as a number.
