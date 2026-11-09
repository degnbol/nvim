# symbol-ls.nvim

Go to code symbols named in markdown code spans. With the cursor in
`` `Foo.bar` ``, `<C-]>` jumps to the definition of `bar` and `K` shows its
hover.

symbol-ls is an in-process LSP server that asks the language servers already
running (lua_ls, basedpyright, rust-analyzer, …) and needs no keymaps of its
own.

Requires Neovim 0.12 and the bundled markdown tree-sitter parsers.

## Setup

```lua
require("symbol-ls").setup()
```

Options, term parsing, ranking and limits: `:help symbol-ls`.

## Tests

```sh
make test
```
