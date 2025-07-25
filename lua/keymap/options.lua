local map = require "utils/keymap"

-- Toggle options.

-- some extra functions that I felt were missing
-- To conceal or not, set by changing the conceallevel between 0
-- and a non-zero value.
local nzConcealLvl
local toggle_conceal = function()
    if vim.opt.conceallevel:get() > 0 then
        nzConcealLvl = vim.opt.conceallevel:get()
        vim.opt.conceallevel = 0
        print("conceallevel = 0")
    else
        vim.opt.conceallevel = nzConcealLvl or 1
        print("conceallevel =", vim.opt.conceallevel:get())
    end
end
local toggle_colcealcursor = function()
    if vim.opt.concealcursor:get():match('n') then
        vim.cmd 'setlocal concealcursor-=n' -- lua version not simple
    else
        vim.opt.concealcursor:append('n')
    end
end
local sidescrolloff
local enable_autoformat = function()
    vim.opt.formatoptions:append('a')
    -- also remove sidescroll offset since there should be enough space on the screen
    -- Keep record of the original value
    sidescrolloff = vim.opt.sidescrolloff:get()
    vim.opt.sidescrolloff = 0
    local notify = "fo+=a | sidescrolloff=0"
    print(notify)
    return notify
end
local disable_autoformat = function()
    vim.opt.formatoptions:remove('a')
    -- reset sidescrolloff.
    -- if the script local var 'sidescrolloff' hasn't been defined
    -- in a call to enable_autoformat, we set it to what is
    -- currently the default in lua/options.lua
    vim.opt.sidescrolloff = sidescrolloff or 12
    local notify = "fo-=a | sidescrolloff=" .. vim.opt.sidescrolloff:get()
    print(notify)
    return notify
end
local autoformat
local enable_wrap = function()
    vim.opt.wrap = true
    -- also disable autoformat when wrapping
    -- Keep record of the original value
    autoformat = vim.opt.formatoptions:get()['a']
    local notify = disable_autoformat() .. " | wrap"
    print(notify)
    return notify
end
local disable_wrap = function()
    vim.opt.wrap = false
    -- reset autoformat.
    -- if the script local var 'autoformat' hasn't been defined
    -- in a call to enable_wrap, we default to false
    local notify
    if autoformat then
        notify = enable_autoformat() .. " | nowrap"
    else
        notify = "nowrap"
    end
    print(notify)
    return notify
end
local toggle_autoformat = function()
    if vim.opt.formatoptions:get()['a'] then
        disable_autoformat()
    else
        enable_autoformat()
    end
end
local function toggle_wrap()
    if vim.opt.wrap:get() then
        disable_wrap()
    else
        enable_wrap()
    end
end
-- %#CursorLineNr# set hl group "CursorLineNr"
-- %## reset hl group, i.e. the group that was active for the line number by default (e.g. set by gitsigns)
-- v:lnum==line(".") line number about to be drawn is cursor line number
-- v:virtnum==0 line number about to be drawn is for buffer line, not virtual line or wrapped line
-- %= right align the following by prepending spaces
-- By setting hl group LineNr on line numbers we overwrite the
-- gitsigns hl group but only for the number itself, which allows
-- us to still use the gitsigns hl group for the following
-- whitespace, which gives the impression there is 1 single line
-- signcolumn.
-- https://old.reddit.com/r/neovim/comments/1dto43b/use_cursorlinenr_highlight_instead_of_gitsigns/
local statuscolumn_relativenumber =
'%#CursorLineNr#%{v:lnum==line(".")&&v:virtnum==0?v:lnum:""}%#LineNr#%=%{v:lnum!=line(".")&&v:virtnum==0?v:relnum:""}%## '
local statuscolumn_number         =
'%#LineNr#%=%{v:lnum!=line(".")&&v:virtnum==0?v:lnum:""}%#CursorLineNr#%{v:lnum==line(".")&&v:virtnum==0?v:lnum:""}%## '
local signcolumn                  = false -- single column signcolumn by itself
local function toggle_signcolumn()
    if signcolumn then
        signcolumn = false
        -- only if presently set by a recent call to toggle_git
        if vim.opt_local.statuscolumn:get() == ' ' then
            vim.opt_local.statuscolumn = ''
        end
        print("no git signcolumn")
    else
        signcolumn = true
        -- only if not already set by number and relativenumber
        if vim.opt_local.statuscolumn:get() == '' then
            vim.opt_local.statuscolumn = ' '
        end
        print("git signcolumn")
    end
end
local function toggle_number()
    if vim.opt_local.number:get() then
        if vim.opt_local.relativenumber:get() then
            vim.opt_local.statuscolumn = statuscolumn_relativenumber
        else
            vim.opt_local.statuscolumn = signcolumn and ' ' or ''
        end
        vim.opt_local.number = false
        print("nonumber")
    else
        if vim.opt_local.relativenumber:get() then
            vim.opt_local.statuscolumn = statuscolumn_relativenumber
        else
            vim.opt_local.statuscolumn = statuscolumn_number
        end
        vim.opt_local.number = true
        print("number")
    end
end
local function toggle_relativenumber()
    if vim.opt_local.relativenumber:get() then
        if vim.opt_local.number:get() then
            vim.opt_local.statuscolumn = statuscolumn_number
        else
            vim.opt_local.statuscolumn = signcolumn and ' ' or ''
        end
        vim.opt_local.relativenumber = false
        print("norelativenumber")
    else
        vim.opt_local.statuscolumn = statuscolumn_relativenumber
        vim.opt_local.relativenumber = true
        print("relativenumber")
    end
end

-- default for yoc is another binding for cursorline which is a lot
-- less useful than conceal
-- TODO: don't toggle cursorline with yo_ and yo-, since it toggles cursor line shown in linenumber, which we always want.
-- Instead make it toggle cursorlineopt +/-= line
map.n('yoc', toggle_conceal, "conceal")
map.n('=sc', toggle_conceal, "conceal")
map.n('yoC', toggle_colcealcursor, "concealcursor")
map.n('=sC', toggle_colcealcursor, "concealcursor")
map.n('yoa', toggle_autoformat, "formatoptions auto")
map.n('=sa', toggle_autoformat, "formatoptions auto")
map.n('<sa', enable_autoformat, "formatoptions auto")
map.n('>sa', disable_autoformat, "formatoptions auto")
map.n('yow', toggle_wrap, "wrap")
map.n('=sw', toggle_wrap, "wrap")
map.n('<sw', enable_wrap, "wrap")
map.n('>sw', disable_wrap, "wrap")
map.n('yon', toggle_number, "number")
map.n('yor', toggle_relativenumber, "relativenumber")
map.n('yog', toggle_signcolumn, "git signcolumn")
