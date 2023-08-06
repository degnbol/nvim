#!/usr/bin/env lua
local rtp = vim.opt.runtimepath:get()[1]
local rtp = rtp .. "/tex/"
dofile(rtp .. "overleaf.lua")
dofile(rtp .. "tables.lua")
dofile(rtp .. "cmds.lua")
dofile(rtp .. "textcolor.lua")

-- also defined in lua/plugins/mini.lua as shift+` but we can just use ` since 
-- I don't think it has use in latex, maybe except in some verbatim code block or something?
require'mini.surround'.config.custom_surroundings['`'] = {
    input = { "``().-()''" },
    output = { left = '``', right = "''" },
}

-- requires kana/vim-textobj-user, see lua/plugins/init.lua
vim.fn["textobj#user#plugin"]("tex", {
    ['latex-ticks'] = {
        pattern = {'``', "''"},
        ['select-a'] = '<buffer> a`',
        ['select-i'] = '<buffer> i`',
    },
})

