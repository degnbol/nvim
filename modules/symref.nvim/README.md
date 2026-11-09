# symref.nvim

Go to code symbols named in markdown code spans. With the cursor in
`` `ShellParse.extract_commands` ``, `<C-]>` jumps to its definition and `K`
shows its hover.

symref is an in-process LSP server that asks the language servers already
running (lua_ls, basedpyright, rust-analyzer, …) and needs no keymaps of its
own.

Requires Neovim 0.12 and the bundled markdown tree-sitter parsers.

## Setup

```lua
require("symref").setup()
```

Options, term parsing, ranking and limits: `:help symref`.

## Tests

```sh
make test
```
