local cmd = vim.cmd

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

-- for guifg, bg

local function fg(group, color)
    cmd("hi " .. group .. " guifg=" .. color)
end

local function bg(group, color)
    cmd("hi " .. group .. " guibg=" .. color)
end

local function fg_bg(group, fgcol, bgcol)
    cmd("hi " .. group .. " guifg=" .. fgcol .. " guibg=" .. bgcol)
end

fg("LineNr", grey)
fg("CursorLineNr", white)
fg("NvimInternalError", red)
fg("VertSplit", line)

-- inactive statuslines as thin splitlines. Not used since statusline is disabled.
-- cmd "hi! StatusLineNC gui=underline guifg=black"
-- cmd "hi StatusLine guibg=black"

-- git signs ---
fg_bg("DiffAdd", green, "none")
fg_bg("DiffAddNr", grey, green)
fg_bg("DiffChange", one_bg2, "none")
fg_bg("DiffChangeNr", grey, nord_blue)
fg_bg("DiffDelete", red, "none")
cmd("hi DiffDeleteNr guifg=" .. grey .. " gui=underline guisp=" .. red)
fg_bg("DiffTopDeleteNr", grey, red)
fg_bg("DiffChangeDeleteNr", grey, dark_purple)
-- file has changed indication
fg_bg("DiffModified", nord_blue, "none")

-- NvimTree
fg("NvimTreeFolderIcon", blue)
fg("NvimTreeFolderName", blue)
fg("NvimTreeIndentMarker", one_bg2)
fg("NvimTreeVertSplit", "black")
-- bg("NvimTreeVertSplit", darker_black)
fg("NvimTreeRootFolder", "black")
bg("NvimTreeNormal", "black")
fg_bg("NvimTreeStatuslineNc", "black", "black")

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
cmd "hi BufferLineBufferVisible guifg=white"
cmd "hi BufferLineNumbersVisible guifg=white"
-- bold instead of italic bold selected file
cmd "hi BufferLineBufferSelected gui=bold"
-- only show tab number for tabs we may want to select, i.e. not current tab
cmd "hi link BufferLineNumbersSelected Ignore"
-- why did this stop working? I had to add this:
cmd "hi BufferLineNumbersSelected guifg=black"


-- dashboard
fg("DashboardHeader", grey_fg)
fg("DashboardCenter", grey_fg)
fg("DashboardShortcut", grey_fg)
fg("DashboardFooter", "black")

bg("MatchParen", "#1f5c6b")

-- Coc indicate error with read undercurl instead of default uncolored underline
cmd("hi CocErrorHighlight cterm=undercurl gui=undercurl guisp=red")

-- HlSearchLens. Less agressive coloring but still matching search hl
cmd "hi HlSearchLens ctermfg=11 ctermbg=242 guifg=#ffdc2d guibg=#5a576e"
cmd "hi HlSearchLensNear ctermfg=11 ctermbg=242 guifg=#19f988 guibg=#5a576e"
