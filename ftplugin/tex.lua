#!/usr/bin/env lua
vim.keymap.set("n", "tsm", "<plug>(vimtex-env-toggle-math)")
vim.keymap.set("n", "dsm", "<plug>(vimtex-env-delete-math)")
vim.keymap.set("n", "csm", "<plug>(vimtex-env-change-math)")
vim.keymap.set("n", "xad", "yaddad", {remap=true})
vim.keymap.set("n", "xid", "yiddid", {remap=true})
-- item with i instead of m and math with m
vim.keymap.set({"o", "x"}, "ai", "<Plug>(vimtex-am)")
vim.keymap.set({"o", "x"}, "ii", "<Plug>(vimtex-im)")
vim.keymap.set({"o", "x"}, "am", "<Plug>(vimtex-a$)")
vim.keymap.set({"o", "x"}, "im", "<Plug>(vimtex-i$)")
-- shorthand to $ just using 4 ($ without shift)
vim.keymap.set({"o", "x"}, "a4", "<Plug>(vimtex-a$)")
vim.keymap.set({"o", "x"}, "i4", "<Plug>(vimtex-i$)")

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

