#!/usr/bin/env lua
local g = vim.g
local utils = require('utils')
local cmd = vim.cmd

-- no default keymappings
g.slime_no_mappings = true
g.slime_target = "kitty"
g.slime_python_ipython = 1 -- we are using ipython https://github.com/jpalardy/vim-slime/tree/main/ftplugin/python


local filetype2command = {
    python="ipython",
    julia="julia",
    r="R"
    -- kitty command will not have access to default setup so doesn't know where R is.
    -- Unfortunately the kitty @ send-text command that is called by slime when using slime_target=kitty sends each line separately where radian then tries to close brackets. 
    -- r="radian --r-binary /Library/Frameworks/R.framework/Resources/R"
}


function kittyWindow()
    -- default to zsh
    ftcommand = filetype2command[vim.bo.filetype] or ""
    fh = io.popen('kitty @ launch --cwd=current --keep-focus ' .. ftcommand)
    window_id = fh:read("*n") -- *n means read number, which means we also strip newline
    fh:close()
    -- set title to the id
    os.execute("kitty @ set-window-title --match id:" .. window_id .. " " .. window_id)
    vim.b.slime_config = {window_id=window_id, listen_on=""}
end

opts = {noremap=true, silent=true}
utils.map("n", "<leader><CR>", ":lua kittyWindow()<CR>", opts)
utils.map("n", "<CR><CR>", ":SlimeSendCurrentLine<CR>j", opts)
-- `> means go to mark named > which will be at the end of the previous selection.
cmd 'xmap <CR> <Plug>SlimeRegionSend()`>'
cmd 'nmap <CR> <Plug>SlimeMotionSend'
-- easily set kitty window id
utils.map("n", "<leader>tk", ':let b:slime_config = {"listen_on": ""}<CR>:let b:slime_config["window_id"] = input("window id: ")<CR>', {noremap=true, silent=true})

