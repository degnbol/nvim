#!/usr/bin/env lua
local luasnip = require "luasnip"

-- https://youtu.be/Dn800rlPIho?t=440
luasnip.config.set_config {
    -- don't jump back into exited snippet
    history = false,
    -- dynamic snippets update as you type
    updateevents = "TextChanged,TextChangedI",
    enable_autosnippets = true,
}
-- load friendly-snippets with luasnip
require("luasnip.loaders.from_vscode").lazy_load()

