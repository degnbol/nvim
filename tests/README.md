# Tests

Unit tests using [plenary.nvim](https://github.com/nvim-lua/plenary.nvim).

**Before writing scratch test files (`/tmp/foo.zsh`, ad-hoc nvim sessions), check `tests/plenary/` for an existing spec covering the area.** New treesitter injections go in `tests/plenary/zsh_injections_spec.lua` (helpers `assert_injection`, `injections_for`); new colorscheme/highlight cases in `colors_spec.lua`; etc. Extend the existing `describe` blocks rather than starting parallel infrastructure.


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

### `ftplugin/` loads, `plugin/` does not

`PlenaryBustedDirectory` (the `make test` target) spawns one nvim per spec file with `--noplugin` (`plenary/test_harness.lua:44,90`), so nothing in `plugin/` runs — a mapping, command or autocmd defined there is simply absent. Source what a spec needs in its describe body: `vim.cmd.runtime("plugin/paths.lua")`.

`PlenaryBustedFile` instead runs in the current nvim, which started without `--noplugin` and so has `plugin/` loaded. A spec that depends on it passes there and fails under `make test`.

`ftplugin/` and the config's filetype detection do apply either way: setting `vim.bo.filetype = "zsh"` gets overridden to `"sh.zsh"`. So setting a filetype can pull in an ftplugin that requires a plugin the harness does not install (`ftplugin/markdown.lua` → `mini.hipatterns`) — `vim.cmd("filetype plugin off")` in the spec avoids that.

### Tree-sitter in tests

`vim.treesitter.get_node()` calls `get_parser(buf)` without an explicit language. This relies on the filetype→parser mapping being correct. Two things to watch:

1. **Compound filetypes**: `"sh.zsh"` maps to the first component `"sh"` by default (no parser installed). The config registers `"zsh"` for `"sh.zsh"` in `plugin/treesitter.lua`, but if that hasn't loaded, add `vim.treesitter.language.register("zsh", "sh.zsh")` in the test file.
2. **Parser must be created before `get_node()`**: Call `vim.treesitter.get_parser(buf):parse()` after setting buffer content. Without this, `get_node()` returns `nil`.
