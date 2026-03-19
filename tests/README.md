# Tests

Unit tests using [plenary.nvim](https://github.com/nvim-lua/plenary.nvim).

## Running tests

Run all tests:
```zsh
nvim --headless -u tests/minimal_init.lua -c "PlenaryBustedDirectory tests/plenary/ {minimal_init = 'tests/minimal_init.lua'}"
```

Run a single test file:
```zsh
nvim --headless -u tests/minimal_init.lua -c "PlenaryBustedFile tests/plenary/colors_spec.lua"
```

## Writing tests

Tests live in `tests/plenary/` and follow the `*_spec.lua` naming convention.

```lua
---@diagnostic disable: undefined-global
describe("module name", function()
    it("does something", function()
        assert.are.equal(expected, actual)
    end)
end)
```

The `---@diagnostic disable: undefined-global` suppresses LSP warnings for busted globals (`describe`, `it`, `assert`, etc.).

## Minimal init

`minimal_init.lua` provides a stripped-down neovim environment for testing:
- Adds the config directory to runtimepath (so `require("utils/...")` works)
- Adds plenary.nvim for the test framework

Tests should avoid triggering the full plugin/ftplugin machinery where possible, as many dependencies won't be loaded.

## Gotchas

### Plugin/ftplugin files still load

Because `minimal_init.lua` adds the config dir to runtimepath, `plugin/` and `ftplugin/` files execute (with errors for missing plugin dependencies, which are non-fatal). This means filetype detection and autocmds from the config are active. For example, setting `vim.bo.filetype = "zsh"` on a buffer will get overridden to `"sh.zsh"` by the config's filetype machinery.

### Tree-sitter in tests

`vim.treesitter.get_node()` calls `get_parser(buf)` without an explicit language. This relies on the filetype→parser mapping being correct. Two things to watch:

1. **Compound filetypes**: `"sh.zsh"` maps to the first component `"sh"` by default (no parser installed). The config registers `"zsh"` for `"sh.zsh"` in `lua/autocmds/treesitter.lua`, but if that hasn't loaded, add `vim.treesitter.language.register("zsh", "sh.zsh")` in the test file.
2. **Parser must be created before `get_node()`**: Call `vim.treesitter.get_parser(buf):parse()` after setting buffer content. Without this, `get_node()` returns `nil`.
