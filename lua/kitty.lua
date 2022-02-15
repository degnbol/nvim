#!/usr/bin/env lua
local g = vim.g
local utils = require'utils'
local cmd = vim.cmd
local json = require'json'
local api = vim.api
local fn = vim.fn

function get_repl()
    return vim.b.repl_id
end
function set_repl(window_id)
    vim.b.repl_id = window_id
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
                                set_repl(win["id"])
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
    -- kitty command will not have access to default setup so doesn't know where R is.
    r="radian --r-binary /Library/Frameworks/R.framework/Resources/R",
    lua="lua",
}

function kittyWindow()
    -- default to zsh
    ftcommand = filetype2command[vim.bo.filetype] or ""
    fh = io.popen('kitty @ launch --cwd=current --keep-focus ' .. ftcommand)
    window_id = fh:read("*n") -- *n means read number, which means we also strip newline
    fh:close()
    -- set title to the id so we can easily set is as target
    os.execute("kitty @ set-window-title --match id:" .. window_id .. " " .. window_id)
    set_repl(window_id)
end


function kittyExists(window_id)
    -- see if kitty window exists by getting text from it and checking if the operation fails (nonzero exit code).
    return os.execute('kitty @ get-text --match id:' .. window_id .. ' > /dev/null 2> /dev/null') == 0
end

function replCheck()
    if get_repl() == nil or not kittyExists(get_repl()) then
        if search_repl() == nil then
            kittyWindow()
        end
    end
end

function kittySend(text)
    replCheck()
    -- pcat.sh uses zsh to do bracketed paste cat from stdin to stdout.
    -- An alternative that fixes indentation but sends each line separately is text:gsub('\n', '\n\x01')
    fh = io.popen('$XDG_CONFIG_HOME/nvim/kittyREPL/kittyPaste.sh ' .. vim.b.repl_id, 'w')
    fh:write(text)
    fh:close()
end

function kittySendLine()
    -- easiest to use stdin rather than putting the text as an arg due to worrying about escaping characters
    kittySend(api.nvim_get_current_line())
end

function kittySendVisual()
    -- gv  = reselect last select (unselected for some reason)
    -- "ky = yank to register k (k for kitty)
    -- `>  = go to mark ">" = end of last visual select
    cmd 'silent normal! gv"ky`>'
    kittySend(fn.getreg('k'))
end

function ReplOperator(type, ...)
    -- `[ = go to start of motion
    -- v or V = select char or lines
    -- `] = go to end of motion
    -- "ky = yank selection to register k
    if type == "char" then
        cmd 'silent normal! `[v`]"ky`]w'
    else -- either type == "line" or "block", the latter I never use
        cmd 'silent normal! `[V`]"ky`]w'
    end
    kittySend(fn.getreg('k'))
end

opts = {noremap=true, silent=true}
utils.map("n", "<leader><CR>", ":lua kittyWindow()<CR>", opts)
utils.map("n", "<CR><CR>", ":lua kittySendLine()<CR>j", opts)
utils.map("x", "<CR>", ":lua kittySendVisual()<CR>", opts)
utils.map('n', "<CR>", 'Operator("v:lua.ReplOperator")', {expr=true, noremap=false})
-- easily set kitty window id for REPL
utils.map("n", "<leader>tt", ':let b:repl_id = input("window id: ")<CR>', {noremap=true, silent=true})

