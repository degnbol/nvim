-- lz.n plugin spec field snippets.
-- Loaded conditionally from ftplugin/lua.lua for lua/plugins/ files.

local ls = require "luasnip"
local s = ls.s
local t = ls.t
local i = ls.i
local c = ls.c
local fmta = require("luasnip.extras.fmt").fmta

return {

-- Trigger fields (lazy-loading)
s({ trig = "event", dscr = "Lazy-load on autocmd event." },
    { t 'event = "', i(1, "DeferredUIEnter"), t '",' }),

s({ trig = "cmd", dscr = "Lazy-load on command." },
    { t 'cmd = ', c(1, {
        fmta('"<>"', { i(1) }),
        fmta('{ "<>" }', { i(1) }),
    }), t ',' }),

s({ trig = "ft", dscr = "Lazy-load on filetype." },
    { t 'ft = ', c(1, {
        fmta('"<>"', { i(1, "lua") }),
        fmta('{ "<>" }', { i(1, "lua") }),
    }), t ',' }),

s({ trig = "keys", dscr = "Lazy-load on key mapping." },
    { t 'keys = ', c(1, {
        fmta([[{
    { "<>", <>, desc = "<>" },
}]], { i(1, "<leader>"), i(2, "function() end"), i(3, "Description") }),
        fmta('"<>"', { i(1, "<leader>") }),
    }), t ',' }),

s({ trig = "colorscheme", dscr = "Lazy-load when colorscheme is set." },
    { t 'colorscheme = "', i(1), t '",' }),

-- Hook functions
s({ trig = "after", dscr = "Runs after the plugin loads." },
    fmta([[after = function()
    <>
end,]], { i(1) })),

s({ trig = "before", dscr = "Runs before the plugin loads." },
    fmta([[before = function()
    <>
end,]], { i(1) })),

s({ trig = "beforeAll", dscr = "Runs before any plugins load (only on first trigger)." },
    fmta([[beforeAll = function()
    <>
end,]], { i(1) })),

-- Base fields
s({ trig = "enabled", dscr = "Conditionally enable/disable the plugin." },
    { t "enabled = ", c(1, { t "false", fmta([[function()
    return <>
end]], { i(1, "true") }) }), t "," }),

s({ trig = "priority", dscr = "Load order for non-lazy plugins. Default: 50, higher = earlier." },
    { t "priority = ", i(1, "1000"), t "," }),

s({ trig = "lazy", dscr = "Force lazy-loading even without triggers." },
    { t "lazy = true," }),

s({ trig = "load", dscr = "Custom load function (overrides global vim.g.lz_n.load)." },
    fmta([[load = function(name)
    <>
end,]], { i(1, "vim.cmd.packadd(name)") })),

}
