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




-- blankline

-- fg("IndentBlanklineChar", line)

-- misc --
fg("LineNr", grey)
fg("NvimInternalError", red)
fg("VertSplit", line)
cmd "hi! EndOfBuffer ctermbg=black ctermfg=black guibg=black guifg=black"

-- inactive statuslines as thin splitlines
cmd "hi! StatusLineNC gui=underline guifg=black"
cmd "hi StatusLine guibg=black"

-- if cursorline is on, only highlight the number not the whole line
cmd "hi clear CursorLine"
fg("cursorlinenr", white)

-- git signs ---
fg_bg("DiffAdd", green, "none")
fg_bg("DiffChange", one_bg2, "none")
fg_bg("DiffChangeNr", nord_blue, "none")
fg_bg("DiffDelete", red, "none")
cmd("hi DiffDeleteNr guifg=" .. grey .. " gui=underline guisp=" .. red)
cmd("hi DiffTopDeleteNr guifg=" .. red)
cmd("hi DiffChangeDeleteNr guifg=" .. dark_purple)
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

fg_bg("BufferLineFill", "black", "black")
fg_bg("BufferLineBackground", "black", "black")
fg_bg("BufferLineBufferVisible", "black", "black")
fg_bg("BufferLineBufferSelected", "black", "black")
-- bold instead of italic bold selected file
cmd "hi BufferLineBufferSelected gui=bold"

-- tabs
fg_bg("BufferLineTab", "black", "black")
fg_bg("BufferLineTabSelected", "black", nord_blue)

fg_bg("BufferLineIndicator", "black", "black")
fg_bg("BufferLineIndicatorSelected", "black", "black")

-- -- separators
fg_bg("BufferLineSeparator", "black", "black")
fg_bg("BufferLineSeparatorVisible", "black", "black")
fg_bg("BufferLineSeparatorSelected", "black", "black")

-- dashboard

fg("DashboardHeader", grey_fg)
fg("DashboardCenter", grey_fg)
fg("DashboardShortcut", grey_fg)
fg("DashboardFooter", "black")

bg("MatchParen", "#1f5c6b")
