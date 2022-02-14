#!/usr/bin/env lua
local g = vim.g
local utils = require'utils'
local cmd = vim.cmd
local json = require'json'

-- no default keymappings
g.slime_no_mappings = true
g.slime_target = "kitty"
g.slime_python_ipython = 1 -- we are using ipython https://github.com/jpalardy/vim-slime/tree/main/ftplugin/python

function slime(window_id)
    vim.b.slime_config = {listen_on="", window_id=window_id}
end

local cmdline2filetype = {
    python3="python",
    ipython="python",
    R="r",
    radian="r",
}

function search_repl()
    fh = io.popen('kitty @ ls')
    json_string = fh:read("*a")
    ls = json.decode(json_string)
    for i_os_win, os_win in ipairs(ls) do
        if os_win["is_focused"] then
            for i_tab, tab in ipairs(os_win["tabs"]) do
                if tab["is_focused"] then
                    for i_win, win in ipairs(tab["windows"]) do
                        cmdline = win["foreground_processes"][1]["cmdline"]
                        -- ["/usr/local/bin/julia", "-t", "4"] -> julia
                        -- [".../R"] -> r
                        -- ["../Python", ".../radian"] -> r
                        -- [".../python3"] -> python
                        for i_arg, arg in ipairs(cmdline) do
                            repl = string.match(arg, "%w+$")
                            repl = cmdline2filetype[repl] or repl
                            if repl == vim.bo.filetype then
                                slime(win["id"])
                                return win["id"]
                            end
                        end
                    end
                end
            end
        end
    end
end

-- register an autocommand to run this when entering buffers
cmd 'au BufEnter * lua search_repl()'

local filetype2command = {
    python="ipython",
    julia="julia",
    r="R",
    -- kitty command will not have access to default setup so doesn't know where R is.
    -- Unfortunately the kitty @ send-text command that is called by slime when using slime_target=kitty sends each line separately where radian then tries to close brackets. 
    -- r="radian --r-binary /Library/Frameworks/R.framework/Resources/R"
    lua="lua",
}

function kittyWindow()
    -- default to zsh
    ftcommand = filetype2command[vim.bo.filetype] or ""
    fh = io.popen('kitty @ launch --cwd=current --keep-focus ' .. ftcommand)
    window_id = fh:read("*n") -- *n means read number, which means we also strip newline
    fh:close()
    -- set title to the id
    os.execute("kitty @ set-window-title --match id:" .. window_id .. " " .. window_id)
    slime(window_id)
end


function slimeCheck()
    if vim.b.slime_config == nil then
        if search_repl() == nil then
            kittyWindow()
        end
    end
end

opts = {noremap=true, silent=true}
utils.map("n", "<leader><CR>", ":lua kittyWindow()<CR>", opts)
utils.map("n", "<CR><CR>", ":lua slimeCheck()<CR>:SlimeSendCurrentLine<CR>j", opts)
-- `> means go to mark named > which will be at the end of the previous selection.
cmd 'xmap <CR> :lua slimeCheck()<CR><Plug>SlimeRegionSend()`>'
cmd 'nmap <CR> <Plug>SlimeMotionSend'
-- easily set kitty window id
utils.map("n", "<leader>tt", ':lua slime(vim.fn.input("window id: "))<CR>', {noremap=true, silent=true})

