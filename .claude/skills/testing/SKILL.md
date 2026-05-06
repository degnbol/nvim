---
name: testing
description: >
  Neovim plugin test isolation and runner patterns. Use when writing or
  debugging plugin tests under modules/, especially with mini.test (busted
  emulation), or when tests pass alone but fail in the full suite.
---

# Testing

Patterns for testing nvim plugins (e.g. `~/nvim/modules/<plugin>.nvim/`).

## mini.test cases re-enter under `vim.wait`

mini.test schedules each case via `vim.schedule(...)`, then runs the queue
on the event loop. `vim.wait(N, condition)` pumps that loop while waiting,
so a queued case can fire **inside** another test's wait. Stubs, autocmds,
module state, and `vim.schedule` callbacks the outer test set up are still
active during the inner case.

Symptom: a test that passes when run alone fails inside the full suite,
with errors that reference state the test never touched (a `vim.uv.fs_stat`
stub from one file leaking into another file's tests is the canonical case).

**Fix: file-level isolation in the runner.** Run each `*.test.lua` in a
fresh `nvim --headless` process. Cross-file pollution becomes structurally
impossible. Within one file tests still share state, but the surface is
small and contained.

```makefile
TEST_FILES := $(shell find lua tests -name "*.test.lua" -o -name "*_test.lua" -o -name "test_*.lua")

test:
	@rc=0; for f in $(TEST_FILES); do \
		nvim --headless -u tests/init.lua -c "lua require('tests.runner').run_file('$$f')" || rc=1; \
	done; \
	exit $$rc
```

Per-file overhead is ~50ms × file count (~2s for 40 files on Apple Silicon).

If a single file still pollutes itself (a `vim.wait`-using helper inside
that file lets the next case re-enter), escalate to **per-test child-neovim**
via `MiniTest.new_child_neovim()`. Spawn a fresh subprocess in `before_each`,
interact via RPC. Slower (~50ms per test) but fully isolates.

## Strip the developer's outer config from rtp

`nvim --headless -u tests/init.lua` does **not** suppress the standard
runtime path. `$XDG_CONFIG_HOME/nvim` (and its `after/`) is still on rtp,
so `plugin/*.lua`, `ftplugin/*.lua`, and autocmds from the developer's
outer nvim config auto-source into the test environment. Tests then
depend on whatever each developer happens to have configured (e.g. a
custom `BufAdd` autocmd reading a missing buffer-local var → spurious
errors during plugin tests that open buffers).

Strip it explicitly in `tests/init.lua`:

```lua
local xdg = os.getenv("XDG_CONFIG_HOME") or (os.getenv("HOME") .. "/.config")
local user_nvim = xdg .. "/nvim"
local user_nvim_after = xdg .. "/nvim/after"
local rtp = vim.api.nvim_get_option_value("runtimepath", {})
local cleaned = {}
for entry in vim.gsplit(rtp, ",", { plain = true }) do
    if entry ~= user_nvim and entry ~= user_nvim_after then
        table.insert(cleaned, entry)
    end
end
vim.api.nvim_set_option_value("runtimepath", table.concat(cleaned, ","), {})
```

Don't overwrite `XDG_CONFIG_HOME` itself in the Makefile — other tools
(`git`'s `core.excludesFile` via `~/.config/git`, etc.) read from it too,
and silently breaking them produces spurious test failures elsewhere.

## Don't depend on the developer's git config

Tests that compare `git ls-files` output (or anything that respects
`--exclude-standard`) inherit the developer's global gitignore via
`core.excludesFile`. A file that's only excluded by *your* global config
will appear in the working set on a clean checkout and break the test.

Fix: commit the project-level `.gitignore` rule explicitly. Don't rely on
"every developer probably has this in their global ignore." Common
offenders: `.claude/settings.local.json`, `.envrc.local`, editor-specific
local files.
