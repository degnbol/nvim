#!/usr/bin/env lua
local g = vim.g
local utils = require('utils')

g.slime_target = "kitty"


local filetype2command = {
    python="ipython",
    julia="julia",
    r="R"
}

function kittyWindow()
    ftcommand = filetype2command[vim.bo.filetype] or ""
    fh = io.popen("kitty @ launch --cwd=current --keep-focus " .. ftcommand)
    window_id = fh:read("*n") -- *n means read number, which means we also strip newline
    fh:close()
    vim.b.slime_config = {window_id=window_id, listen_on=""}
end

utils.map("n", "<leader>k", ":lua kittyWindow()<CR>", {noremap = true, silent = true})


