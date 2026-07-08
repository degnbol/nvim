---
name: plugins
description: Plugin management in this config — vim.pack (install/update/lockfile) + lz.n (lazy-loading). Use when editing lua/pack_specs.lua, lua/pack_hooks.lua, lua/plugins/*.lua, adding/removing/disabling plugins, or debugging plugin load order and startup overhead.
---

# Package Management: vim.pack + lz.n

## Architecture

Plugins are installed by nvim 0.12's built-in `vim.pack` and lazy-loaded by `lz.n`.

- **`vim.pack.add(specs)`** — clones/updates plugins into `~/.local/share/nvim/site/pack/core/opt/`
- **`lz.n`** — registers triggers (event/cmd/ft/keys) and calls `vim.cmd.packadd(name)` on demand

## vim.pack.add behaviour during init.lua

`vim.pack.add()` calls `:packadd!` (with bang) for every plugin during init.lua sourcing.
`:packadd!` adds the plugin directory to `rtp` but does NOT source `plugin/` files.

However, neovim's startup sequence has a "loading rtp plugins" phase after init.lua
completes, which sources ALL `plugin/*.vim` and `plugin/*.lua` files found on rtp. This
means **every plugin registered with vim.pack.add ends up on rtp and has its plugin/
files sourced at startup**, regardless of lz.n lazy-loading triggers.

### Consequence for lz.n

lz.n's lazy loading only controls when `packadd` (without bang) runs and when
`before`/`after` hooks execute. But since vim.pack already added all plugins to rtp
via `:packadd!`, the plugins' own `plugin/` directory files run during the rtp scan
phase anyway. lz.n effectively only defers the `after` hook (setup/config), not the
plugin's own initialisation files.

This is a fundamental mismatch between vim.pack (which puts everything on rtp) and
lz.n (which assumes plugins are NOT on rtp until explicitly loaded).

### The `load` option

`vim.pack.add(specs, opts)` accepts a `load` option:

| Value | Behaviour |
|-------|-----------|
| `false` (default during init.lua) | `:packadd!` — adds to rtp, no plugin/ sourcing |
| `true` (default after init) | `:packadd` — adds to rtp AND sources plugin/ files |
| `function(plug_data)` | Custom loader — vim.pack does NOT call packadd at all |

Using `load = function() end` prevents vim.pack from touching rtp entirely, letting
lz.n be the sole mechanism for loading plugins. Plugins without lz.n specs then need
explicit `packadd` calls or lz.n specs to ever load.

### Diagnosing startup overhead

Use `--startuptime` on the **Embedded** process (second section in output):

```zsh
nvim --startuptime timing.txt --headless +qa
```

Check how many opt plugins end up on rtp:

```lua
nvim --headless +'lua local rtp = vim.opt.rtp:get(); local n = 0; for _, p in ipairs(rtp) do if p:find("opt/") then n = n + 1 end end; print(n)' +qa
```

## lz.n spec fields

| Field | Effect |
|-------|--------|
| `event` | Load on nvim event (e.g. `"InsertEnter"`, `"DeferredUIEnter"`) |
| `cmd` | Load on command (e.g. `"Telescope"`) |
| `ft` | Load on filetype |
| `keys` | Load on keypress |
| `lazy = true` | Don't load eagerly (redundant if any trigger is set) |
| `enabled = false` | Skip entirely (no packadd, no hooks) |
| `before` | Runs before packadd (keymaps, require wrappers) |
| `after` | Runs after packadd (setup/config calls) |
| `load` | Custom load function (replaces default packadd) |
| `priority` | Numeric, higher loads first among eager plugins (default 50) |

### How lz.n determines laziness

A plugin is lazy if `spec.lazy == true` OR any handler field (event/cmd/ft/keys/colorscheme) is present. Plugins without triggers and without `lazy = true` load eagerly during `lz.n.load()`.

### DeferredUIEnter

Fires after UIEnter via `vim.schedule`. Equivalent to lazy.nvim's `VeryLazy`. Still
runs during `--startuptime` measurement — it defers past the UI attach but the user
still waits for it.

### trigger_load gotcha

`trigger_load("name")` only finds plugins registered with a handler. Plugins with
`lazy = true` but no trigger field (event/cmd/ft/keys) are invisible to trigger_load.
Use `vim.cmd.packadd("name")` + setup in the parent's hooks instead.

## Eager keymaps for lazy plugins

When plugin keymaps must exist at startup (e.g. for which-key/mini.clue discovery)
but the plugin itself should lazy-load, register the keymaps outside the lz.n spec
(at module top-level in the plugins/*.lua file) with callbacks that trigger loading:

```lua
-- At module level (runs when lz.n discovers the file):
local function load_picker()
    require("lz.n").trigger_load("fzf-lua")
end

map.n("<leader>fg", function()
    load_picker()                          -- packadd + before/after hooks
    require("fzf-lua").grep()              -- now safe to require
end, "Grep")

return {
    {
        "fzf-lua",
        lazy = true,
        cmd = "FzfLua",                    -- handler needed for trigger_load
        after = function() ... end,
    },
}
```

Key points:
- **Keymaps in `before`/`after` only run when the plugin loads** — they won't exist
  at startup, so which-key menus won't show them and prefix keys (e.g. `<leader>f`)
  may time out or conflict.
- **Module-level code in plugins/*.lua runs eagerly** when lz.n requires the file to
  read specs. Use this for keymap registration.
- **`require("plugin")` fails for opt packages** that haven't been `packadd`'d yet.
  Always call `trigger_load()` or `vim.cmd.packadd()` before the first `require()`.
- **`trigger_load` needs a handler** (cmd/event/ft/keys) on the spec. Plugins with
  only `lazy = true` are invisible to it — use `vim.cmd.packadd()` directly instead.
- **`trigger_load` is idempotent** — calling it after the plugin is already loaded is
  a no-op, so there's no cost to calling it in every keymap callback.

## `after` runs once per plugin load, not per buffer

`before` and `after` execute a single time — when lz.n first triggers packadd for the
plugin. They do NOT re-run for each new buffer that matches the `ft` trigger or each
keypress that matches the `keys` trigger. Placing buffer-local setup inside `after`
silently breaks on the second-and-later buffer:

```lua
-- BUG: only the first asciidoc buffer gets these.
{
    "vim-asciidoc",
    ft = "asciidoc",
    after = function()
        vim.cmd [[compiler asciidoctor2pdf]]                          -- window-local
        vim.keymap.set("n", "<leader>cc", ..., { buffer = true })     -- buffer-local
        vim.opt_local.comments = "://"                                -- buffer-local
        vim.api.nvim_create_autocmd("ColorScheme", { buffer = 0, ... }) -- pinned to buf
    end,
}
```

**Fix**: keep only genuinely global setup in `after` (`vim.g.*`, plugin `setup()`,
global `ColorScheme` autocmds via `hi.onColorScheme` or equivalent). Move per-buffer
settings to `ftplugin/<ft>.lua` or `after/ftplugin/<ft>.lua`, which neovim's
`filetypeplugin` re-runs for every matching buffer. Use `after/ftplugin/` when the
settings depend on the lz.n-loaded plugin's own ftplugin having run first (e.g. to
override a `comments=` clear, or call `compiler <name>` once the plugin's
`compiler/<name>.vim` is available).

## File map

- `lua/pack_specs.lua` — `vim.pack.add()` registry: declares *what* to install (all remote plugins)
- `lua/pack_hooks.lua` — `PackChanged` build hooks (mason, treesitter, etc.)
- `lua/plugins/*.lua` — lz.n specs: declares *how/when* to load; auto-discovered via `require("lz.n").load("plugins")`
- **Dev plugins** (`modules/`) — added to rtp manually in `init.lua`; their lz.n specs use `load = function() end` so vim.pack doesn't touch rtp for them

Plugins in `pack_specs.lua` without a corresponding lz.n spec are installed to opt/ but
have no loading mechanism. With the current vim.pack behaviour, they still end up on
rtp and their plugin/ files still run. With a custom `load` function, they would
genuinely never load unless explicitly packadd'd.

## Disabling plugins

Use `enabled = false` in the lz.n spec. Keep the full spec (setup config, keymaps,
etc.) intact so the plugin can be re-enabled later by removing the flag. Never delete
the spec body when disabling.
