#!/usr/bin/env lua

-- is lazy.nvim plugin
if vim.api.nvim_buf_get_name(0):match("lua/plugins/") then
    local ls = require("luasnip")
    -- https://github.com/L3MON4D3/LuaSnip/blob/master/DOC.md
    local s = ls.s
    local t = ls.t
    local i = ls.i
    local c = ls.c
    local fmta = require("luasnip.extras.fmt").fmta
    
    ls.add_snippets("lua", {
-- https://github.com/folke/lazy.nvim

s({trig="dir", dscr="Path to local plugin."},
{t'dir = "$XDG_CONFIG_HOME/nvim/', i(1, "NAME"), t'.nvim",'}),

s({trig="url", dscr="A custom git url where the plugin is hosted."},
{t'url = "https://', i(1, "gitlab.com/USER/REPO"), t'",'}),

s({trig="name", dscr="A custom name for the plugin used for the local plugin directory and as the display name."},
{t'name = "', i(1, "NAME"), t'",'}),

s({trig="dev", dscr="When true, a local plugin directory will be used instead. See `config.dev`."},
{t"dev = true,"}),

s({trig="lazy", dscr="When true, the plugin will only be loaded when needed. Lazy-loaded plugins are automatically loaded when their Lua modules are required, or when one of the lazy-loading handlers triggers."},
{t"lazy = true,"}),

s({trig="enabled", dscr="When false, or if the function returns false, then this plugin will not be included in the spec."},
{t"enabled = false,"}),

s({trig="cond", dscr="When false, or if the function returns false, then this plugin will not be loaded. Useful to disable some plugins in vscode, or firenvim for example."},
{t"cond = ", i(1, "false"), t","}),

s({trig="dependencies", dscr="A list of plugin names or plugin specs that should be loaded when the plugin loads. Dependencies are always lazy-loaded unless specified otherwise. When specifying a name, make sure the plugin spec has been defined somewhere else."},
{t'dependencies = { "', i(1, "USER/REPO"), t'" },'}),

s({trig="init", dscr="`init` functions are always executed during startup."},
fmta([[init = function () {
    <>
},]], {i(1)})),

s({trig="opts", dscr="`opts` should be a table (will be merged with parent specs), return a table (replaces parent specs) or should change a table. The table will be passed to the `Plugin.config()` function. Setting this value will imply `Plugin.config()`."},
fmta([[opts = {
    <>
},]], {i(1)})),

s({trig="config", dscr="`config` is executed when the plugin loads. The default implementation will automatically run `require(MAIN).setup(opts)`. Lazy uses several heuristics to determine the plugin's `MAIN` module automatically based on the plugin's name. See also `opts`. To use the default implementation without `opts` set config to `true`."},
{t"config = ", c(1, {fmta([[function ()
    <>
end]], {i(1, 'require"PLUGIN".setup {}')}), t"true"}), t","}),

s({trig="main", dscr="You can specify the `main` module to use for `config()` and `opts()`, in case it can not be determined automatically. See `config()`."},
{t'main = "', i(1, "MAIN"), t'",'}),

s({trig="build", dscr="`build` is executed when a plugin is installed or updated. Before running `build`, a plugin is first loaded. If it's a string it will be ran as a shell command. When prefixed with `:` it is a Neovim command. You can also specify a list to executed multiple build commands. Some plugins provide their own `build.lua` which is automatically used by lazy. So no need to specify a build step for those plugins."},
{t'build = "', i(1, "make"), t'",'}),

s({trig="branch", dscr="Branch of the repository."},
{t'branch = "', i(1, "dev"), t'",'}),

s({trig="tag", dscr="Tag of the repository."},
{t'tag = "', i(1, "legacy"), t'",'}),

s({trig="commit", dscr="Commit of the repository."},
{t'commit = "', i(1), t'",'}),

s({trig="version", dscr="Version to use from the repository. Full `Semver <https://devhints.io/semver>` ranges are supported."},
{t"version = false,"}),

s({trig="pin", dscr="When `true`, this plugin will not be included in updates."},
{t"pin = true,"}),

s({trig="submodules", dscr="When `false`, git submodules will not be fetched. Defaults to `true`."},
{t"submodules = false,"}),

s({trig="event", dscr="Lazy-load on event. Events can be specified as `BufEnter` or with a pattern like `BufEnter *.lua`."},
{t'event = "', i(1), t'",'}),

s({trig="cmd", dscr="Lazy-load on command."},
{t'cmd = "', i(1), t'",'}),

s({trig="ft", dscr="Lazy-load on filetype."},
{t'ft = "', i(1, "lua"), t'",'}),

s({trig="keys", dscr="Lazy-load on key mapping."},
{t'keys = "', i(1), t'",'}),

s({trig="module", dscr="Do not automatically load this module when it's required somewhere."},
{t"module = false,"}),

s({trig="priority", dscr="Only useful for start plugins (`lazy=false`) to force loading certain plugins first. Default priority is `50`. It's recommended to set this to a high number for colorschemes."},
{t"priority = ", i(1, "50"), t","}),

s({trig="optional", dscr="When a spec is tagged optional, it will only be included in the final spec, when the same plugin has been specified at least once somewhere else without `optional`. This is mainly useful for Neovim distros, to allow setting options on plugins that may/may not be part of the user's plugins."},
{t"optional = true,"}),

    })
end
