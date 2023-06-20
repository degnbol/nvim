#!/usr/bin/env lua

local function hi(name, val)
    vim.api.nvim_set_hl(0, name, val)
end
local function fg(name, color)
    hi(name, {fg=color})
end
local function bg(name, color)
    hi(name, {bg=color})
end
local function fgbg(name, fgcol, bgcol)
    hi(name, {fg=fgcol, bg=bgcol})
end
local function link(name, linkto)
    hi(name, {link=linkto})
end
local function rev(name)
    hi(name, {reverse=true})
end

-- Italic highlight group doesn't actually make terminal text italic by default.
hi("Italic", {italic=true})

-- set completion menu bg to main bg and make scrollbar minimal
bg("Pmenu", nil)
bg("PmenuSbar", nil)
rev("PmenuSel")

-- if you want, you can also show search as flipping foreground and background
-- note, this doesn't break the green IncSearch
rev("Search")

-- @variable seemed to start being linked to @identifier after update?
hi("@variable", {gui=nil})

-- remove bg on concealchars
bg("Conceal", nil)


---- bufferline
-- bold instead of italic+bold selected
hi("BufferLineBufferSelected", {bold=true})

-- hide weird error highlights on cmp menu popup for docs for functions.
-- When checking the buftype (normal KK:echo &buftype) it says "nofile" so I 
-- guess that is set after running this file.
-- This means, if you have any problems with lacking error highlight in a real 
-- file it means you have to replace this if with an aucmd that is called when 
-- a buffer is set as nofile.
if vim.bo.buftype == "" and vim.bo.filetype == "" then
    link("Error", "Ignore")
end

