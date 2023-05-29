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

local colors = require "themes/onedark"
local white = colors.white
local darker_black = colors.darker_black
local black = colors.black
local black2 = colors.black2
local one_bg = colors.one_bg
local one_bg2 = colors.one_bg2
local one_bg3 = colors.one_bg3
local light_grey = colors.light_grey
local grey = colors.grey
local grey_fg = colors.grey_fg
local red = colors.red
local line = colors.line
local green = colors.green
local nord_blue = colors.nord_blue
local blue = colors.blue
local yellow = colors.yellow
local purple = colors.purple
local dark_purple = colors.dark_purple

-- Italic highlight group doesn't actually make terminal text italic by default.
hi("Italic", {italic=true})

fgbg("LineNr", grey, "black")
fgbg("CursorLineNr", white, "black")
fg("NvimInternalError", red)
fgbg("VertSplit", line, "black")

-- inactive statuslines as thin splitlines. Not used since statusline is disabled.
-- hi("StatusLineNC", {underline=true, fg="black"})
-- bg("StatusLine", "black")

-- git signs ---
fgbg("DiffAdd", green, nil)
fgbg("DiffAddNr", grey, green)
fgbg("DiffChange", one_bg2, nil)
fgbg("DiffChangeNr", grey, nord_blue)
fgbg("DiffDelete", red, nil)
hi("DiffDeleteNr", {fg=grey, underline=true, special=red})
fgbg("DiffTopDeleteNr", grey, red)
fgbg("DiffChangeDeleteNr", grey, dark_purple)
-- file has changed indication
fgbg("DiffModified", nord_blue, nil)

-- NvimTree
fg("NvimTreeFolderIcon", blue)
fg("NvimTreeFolderName", blue)
fg("NvimTreeIndentMarker", one_bg2)
fg("NvimTreeVertSplit", "black")
-- bg("NvimTreeVertSplit", darker_black)
fg("NvimTreeRootFolder", "black")
bg("NvimTreeNormal", "black")
fgbg("NvimTreeStatuslineNc", "black", "black")

-- telescope
fg("TelescopeBorder", grey)
fg("TelescopePromptBorder", grey)
fg("TelescopeResultsBorder", grey)
fg("TelescopePreviewBorder", grey)

-- LspDiagnostics ---

-- error / warnings
fg("LspDiagnosticsSignError", red)
fg("LspDiagnosticsVirtualTextError", red)
fg("LspDiagnosticsSignWarning", yellow)
fg("LspDiagnosticsVirtualTextWarning", yellow)

-- info
fg("LspDiagnosticsSignInformation", green)
fg("LspDiagnosticsVirtualTextInformation", green)

-- hint
fg("LspDiagnosticsSignHint", purple)
fg("LspDiagnosticsVirtualTextHint", purple)

-- bufferline
fg("BufferLineBackground", colors.grey_fg2)
fg("BufferLineNumbers", colors.grey_fg2)
fg("BufferLineBufferVisible", "white")
fg("BufferLineNumbersVisible", "white")
-- bold instead of italic bold selected file
hi("BufferLineBufferSelected", {bold=true})
-- only show tab number for tabs we may want to select, i.e. not current tab
link("BufferLineNumbersSelected", "Ignore")
-- why did this stop working? I had to add this:
fg("BufferLineNumbersSelected", "black")

-- dashboard
fg("DashboardHeader", grey_fg)
fg("DashboardCenter", grey_fg)
fg("DashboardShortcut", grey_fg)
fg("DashboardFooter", "black")

bg("MatchParen", "#1f5c6b")

-- HlSearchLens. Less agressive coloring but still matching search hl
hi("HlSearchLens", {ctermfg=11, ctermbg=242, fg="#ffdc2d", bg="#5a576e"})
hi("HlSearchLensNear", {ctermfg=11, ctermbg=242, fg="#19f988", bg="#5a576e"})

-- set completion menu bg to main bg and make scrollbar minimal
link("Pmenu", "Normal")
link("PmenuSbar", "Ignore")
link("PmenuThumb", "Visual")

bg("Folded", "black")
bg("FoldColumn", "black")

hi("Operator", {ctermfg=5, fg="#ae94f9"})

fg("NonText", light_grey)

-- if you want, you can also show search as flipping foreground and background
-- note, this doesn't break the green IncSearch
hi("Search", {reverse=true, fg=nil, bg=nil})

-- @variable seemed to start being linked to @identifier after update?
hi("@variable", {gui=nil})

-- dim the source of completion item.
-- This column won't be visible if they are all "", see lua/plugins/completion.lua
fg("CmpItemMenu", light_grey)
