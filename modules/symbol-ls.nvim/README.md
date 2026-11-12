# symbol-ls.nvim

Go to code symbols named in markdown code spans. With the cursor in
`` `Foo.bar` ``, `<C-]>` jumps to the definition of `bar` and `K` shows its
hover. In Visual mode, `K` and `<C-]>` look up the selected text instead.

symbol-ls is an in-process LSP server that asks the language servers already
running (lua_ls, basedpyright, rust-analyzer, …), else runs tree-sitter `tags`
queries on the git-tracked files that contain the symbol. Normal mode uses the
core LSP keymaps. symbol-ls maps only Visual `K` and `<C-]>`.

Requires Neovim 0.12 and the bundled markdown tree-sitter parsers. The
`tags` fallback also needs git, curl and nvim-treesitter.

## Setup

```lua
require("symbol-ls").setup()
```

Options, term parsing, ranking and limits: `:help symbol-ls`.

## Tests

```sh
make test
```
