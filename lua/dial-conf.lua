#!/usr/bin/env lua
-- https://github.com/monaqa/dial.nvim

local augend = require("dial.augend")
require("dial.config").augends:register_group {
    default = {
        augend.integer.alias.decimal_int,
        augend.constant.alias.bool,    -- boolean value (true <-> false)
        augend.constant.new{ elements={"True", "False"}, word=true, cyclic=true, }, -- python
        augend.hexcolor.new{ case="lower", },
        augend.date.alias["%Y/%m/%d"],
        augend.date.alias["%Y-%m-%d"],
        augend.date.alias["%m/%d"],
        augend.date.alias["%H:%M"],
    },
}


vim.api.nvim_set_keymap("n", "<C-a>", require("dial.map").inc_normal(), {noremap = true})
vim.api.nvim_set_keymap("n", "<C-x>", require("dial.map").dec_normal(), {noremap = true})
vim.api.nvim_set_keymap("v", "<C-a>", require("dial.map").inc_visual(), {noremap = true})
vim.api.nvim_set_keymap("v", "<C-x>", require("dial.map").dec_visual(), {noremap = true})
vim.api.nvim_set_keymap("v", "g<C-a>", require("dial.map").inc_gvisual(), {noremap = true})
vim.api.nvim_set_keymap("v", "g<C-x>", require("dial.map").dec_gvisual(), {noremap = true})

