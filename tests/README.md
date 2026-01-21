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
