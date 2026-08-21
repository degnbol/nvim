# kittyPath: dangling syntax reference

`after/syntax/kitty.vim` ends its only line with `nextgroup=kittyPath`, but no
`kittyPath` syntax item exists — verified with `syntax list kittyPath` (empty),
`nvim_get_hl(0, {name = "kittyPath"})` (`{}`), and `synstack()` over a
`startup_session` path, which reports `kittySt` rather than `kittyPath`. The
`nextgroup=` half of that line has therefore never fired under vim.pack.

The group used to be defined by an *uncommitted in-place edit* to the plugin's
own `syntax/kitty.vim` in the lazy.nvim checkout. The vim.pack copy is a clean
checkout, so the migration dropped it. As of writing, the only copy is that
uncommitted edit under `~/.local/share/nvim/lazy/vim-kitty/`, which is slated
for deletion — hence recorded here verbatim:

```diff
-syn match kittyInclude '^\s*\(env\|glob\)\?include' display
+syn match kittyPath '.*$' contained
+syn match kittyInclude '^\s*\(glob\|gen\)\?include' display nextgroup=kittyPath
+syn match kittyInclude '^\s*envinclude' display

 hi def link kittyNumber Number
+hi def link kittyPath String
```

Two halves worth separating: the `kittyPath` group itself, and the split of
`envinclude` out of the `\(env\|glob\)\?include` alternation so that only
`include`/`globinclude`/`geninclude` take a path argument. Both `envinclude`
and `geninclude` are real keywords — see
`modules/kitty-conf.nvim/lua/kitty-conf/kitty_options.json`.

**Fix belongs in `after/syntax/kitty.vim`, not the plugin.** Editing the
plugin's own file is what lost it: any plugin update overwrites it. Overriding
from `after/` needs `syn clear kittyInclude` before redefining, since the
plugin's rule has already been read by then (compound-filetype load order:
neovim skill → "Compound Filetypes and Runtime Loading").

Left undone deliberately — a mis-ported syntax rule silently degrades
highlighting, so this wants a visual check against a real `kitty.conf`.
