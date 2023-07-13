#!/usr/bin/env lua
-- some sane defaults that a colorscheme may overwrite

local function hl(name, val)
    vim.api.nvim_set_hl(0, name, val)
end
local function fg(name, color)
    hl(name, {fg=color})
end
local function bg(name, color)
    hl(name, {bg=color})
end
local function fgbg(name, fgcol, bgcol)
    hl(name, {fg=fgcol, bg=bgcol})
end
local function link(name, linkto)
    hl(name, {link=linkto})
end
local function rev(name)
    hl(name, {reverse=true})
end

function gethl(name)
    return vim.api.nvim_get_hl(0, {name=name, link=false})
end
-- update subset of settings for a highlight group instead of replacing them all
function modhl(name, val)
    local hls = gethl(name)
    hls = vim.tbl_extend("force", hls, val)
    hl(name, hls)
end

-- Italic highlight group doesn't actually make terminal text italic by default.
hl("Italic", {italic=true})

-- remove bg on concealchars
bg("Conceal", nil)

-- set completion menu bg to main bg and make scrollbar minimal
bg("Pmenu", nil)
bg("PmenuSbar", nil)
rev("PmenuSel")

fg("LineNr", "grey")
fg("CursorLineNr", nil)

-- hide weird error highlights on cmp menu popup for docs for functions.
-- When checking the buftype (normal KK:echo &buftype) it says "nofile" so I 
-- guess that is set after running this file.
-- This means, if you have any problems with lacking error highlight in a real 
-- file it means you have to replace this if with an aucmd that is called when 
-- a buffer is set as nofile.
if vim.bo.buftype == "" and vim.bo.filetype == "" then
    link("Error", "Ignore")
end

-- @variable seemed to start being linked to @identifier after update?
hl("@variable", {gui=nil})

---- git signs link to Gitsigns equivalent which gets defined
link("DiffAddNr", "GitsignsAddLn")
link("DiffChangeNr", "GitsignsChangeLn")
-- DiffDeleteNr is set in afterColorscheme since it needs to know a 
-- color defined by Gitsigns and can't just link
link("DiffTopDeleteNr", "GitsignsDeleteVirtLn") -- edge-case where first line(s) of file is deleted
link("DiffChangeDeleteNr", "GitsignsChangedeleteLn")

---- bufferline
-- bold instead of italic+bold selected
hl("BufferLineBufferSelected", {bold=true})


local function afterColorscheme()
    -- fix issue where colorscheme change removes all telescope highlight groups
    vim.cmd 'silent Lazy reload telescope.nvim'

    -- special decides the color for the underline
    hl("DiffDeleteNr", {underline=true, special=gethl("DiffDelete")["fg"], fg=nil})

    local spellBad = gethl("spellBad")
    -- the default already has red undercurl so we don't wanna mess with that but 
    -- other themes sets the fg instead
    if spellBad["special"] then
        spellBad = spellBad["special"]
    else
        spellBad = spellBad["fg"]
    end
    hl("SpellBad", {undercurl=true, special=spellBad})

    -- replace whatever search was doing completely but only modify IncSearch's 
    -- reverse setting so it remains a different color from Search, probably 
    -- green.
    rev("Search")
    modhl("IncSearch", {reverse=true})

    -- never italic comments but italize builtin stuff
    -- :h group-name
    -- Italics are not default since terms aren't assumed to support it,
    -- however bold is sometimes set for highlight groups for some languages.
    -- italic for builtin reserved words makes sense since they are few and 
    -- don't need to be too legible.
    -- highlight groups starting with @ are from treesitter and others from regex 
    -- syntax groups. Treesitter highlights take priority over regex.
    -- Which one is edited matters, since editing one with color 
    -- will combine the color with the italic, while editing the other may replace 
    -- the color with italic. Before making changes look at test files in 
    -- testfiles/
    modhl("Comment", {italic=false})
    modhl("Operator", {bold=true})
    modhl("Include", {italic=true})
    modhl("Repeat", {italic=true})
    modhl("Label", {italic=true})
    modhl("Conditional", {italic=true})
    modhl("Exception", {italic=true})
    modhl("@include", {italic=true})
    modhl("Keyword", {italic=true, bold=false})
    modhl("@keyword", {italic=true, bold=false})
    modhl("@keyword.function", {italic=true})
    modhl("@keyword.return", {italic=true})
    modhl("@keyword.operator", {italic=true})
    modhl("@parameter", {italic=false})

    -- NonText shouldn't be exactly like comments
    if gethl("Comment")['fg'] == gethl("NonText")['fg'] then
        fg("NonText", "gray")
    end
end

local grp = vim.api.nvim_create_augroup("afterColorscheme", {clear=true})
vim.api.nvim_create_autocmd("Colorscheme", {
    pattern = "*", group = grp,
    callback = afterColorscheme,
})
-- Not sure why it needs a 0 ms delay
vim.api.nvim_create_autocmd("VimEnter", {
    pattern = "*", group = grp,
    callback = function () vim.defer_fn(afterColorscheme, 0) end
})

