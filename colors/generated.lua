-- Generated neovim colorscheme from terafox-purple
-- Do not edit — regenerate with: ./generate.sh themes/terafox-purple.toml

vim.o.background = "dark"
vim.g.colors_name = "generated"

local hi = vim.api.nvim_set_hl

-- UI
hi(0, "Normal", { bg = "#152528", fg = "#e6eaea" })
hi(0, "NormalFloat", { bg = "#101718", fg = "#e6eaea" })
hi(0, "NormalNC", { bg = "#152528", fg = "#e6eaea" })
hi(0, "Visual", { bg = "#293e40" })
hi(0, "CursorLine", { bg = "#273D41" })
hi(0, "CursorLineNr", { bold = true, fg = "#fda47f" })
hi(0, "LineNr", { fg = "#537479" })
hi(0, "SignColumn", { fg = "#537479" })
hi(0, "ColorColumn", { bg = "#1D3033" })
hi(0, "Pmenu", { bg = "#293e40", fg = "#e6eaea" })
hi(0, "PmenuSel", { bg = "#4C696B" })
hi(0, "PmenuSbar", { bg = "#293e40", fg = "#e6eaea" })
hi(0, "PmenuThumb", { bg = "#4C696B" })
hi(0, "FloatBorder", { fg = "#537479" })
hi(0, "WinSeparator", { fg = "#000000" })
hi(0, "VertSplit", { fg = "#000000" })
hi(0, "StatusLine", { fg = "#101718" })
hi(0, "StatusLineNC", { fg = "#101718" })
hi(0, "TabLine", { bg = "#1D3033", fg = "gray" })
hi(0, "TabLineFill", { bg = "#101718", fg = "gray" })
hi(0, "TabLineSel", { bg = "#537479", fg = "#152528" })
hi(0, "Folded", { bold = true, fg = "#304A4F" })
hi(0, "MoreMsg", { link = "NonText" })
hi(0, "Question", { link = "NonText" })
hi(0, "MsgArea", { bg = "#101718" })
hi(0, "NonText", { bold = true, fg = "#304A4F" })
hi(0, "MatchParen", { underline = true })
hi(0, "IncSearch", { bold = true, standout = true })
hi(0, "Search", { fg = "gray", standout = true })
hi(0, "CurSearch", { link = "Search" })
hi(0, "Directory", { bold = true })
hi(0, "Cursor", { bg = "#e6eaea", fg = "#152528" })
hi(0, "TermCursor", { bg = "#e6eaea", fg = "#152528" })

-- Diff
hi(0, "DiffAdd", { bg = "#294145" })
hi(0, "DiffDelete", { bg = "#29383A", fg = "#152528" })
hi(0, "DiffChange", { bg = "#253C41" })
hi(0, "DiffText", { bg = "#2F4C52" })

-- Syntax
hi(0, "Comment", { fg = "#4e5157", italic = false })
hi(0, "String", { fg = "#7aa4a1" })
hi(0, "Function", { fg = "#5a93aa" })
hi(0, "Keyword", { fg = "#8577af", italic = true })
hi(0, "Statement", { fg = "#8577af" })
hi(0, "Conditional", { fg = "#8577af", italic = true })
hi(0, "Repeat", { fg = "#8577af", italic = true })
hi(0, "Exception", { fg = "#8577af", italic = true })
hi(0, "Label", { bold = true, fg = "#8577af", italic = false })
hi(0, "Operator", { bold = true, fg = "#5a93aa" })
hi(0, "Type", { fg = "#a1cdd8" })
hi(0, "Number", { fg = "#7B9FA2" })
hi(0, "Boolean", { fg = "#7B9FA2", italic = true })
hi(0, "Constant", { fg = "#7B9FA2" })
hi(0, "Identifier", { fg = "#e6eaea" })
hi(0, "PreProc", { fg = "#8577af" })
hi(0, "Include", { link = "PreProc" })
hi(0, "Special", { fg = "#5a93aa" })
hi(0, "Macro", { fg = "#8eb2af" })
hi(0, "Delimiter", { fg = "#8577af" })
hi(0, "Title", { bold = true, underdouble = true })
hi(0, "Bold", { bold = true })
hi(0, "Italic", { italic = true })
hi(0, "Underlined", { underline = true })

-- Diagnostics
hi(0, "DiagnosticError", { fg = "#eb746b" })
hi(0, "DiagnosticWarn", { fg = "#fda47f" })
hi(0, "DiagnosticInfo", { fg = "#5a93aa" })
hi(0, "DiagnosticHint", { fg = "#7aa4a1" })
hi(0, "DiagnosticUnderlineError", { sp = "#eb746b", undercurl = true })
hi(0, "DiagnosticUnderlineWarn", { sp = "#fda47f", underdotted = true })
hi(0, "DiagnosticUnderlineInfo", { sp = "#5a93aa", underdotted = true })
hi(0, "DiagnosticUnderlineHint", { sp = "#7aa4a1", underdotted = true })
hi(0, "ErrorMsg", { fg = "#eb746b" })
hi(0, "WarningMsg", { fg = "#fda47f" })
hi(0, "Error", { sp = "#eb746b", undercurl = true })
hi(0, "SpellBad", { sp = "#eb746b", undercurl = true })

-- Treesitter
hi(0, "@comment", { link = "Comment" })
hi(0, "@string", { link = "String" })
hi(0, "@string.escape", { fg = "#8eb2af" })
hi(0, "@string.documentation", { fg = "#636B72" })
hi(0, "@string.regexp", { link = "@string" })
hi(0, "@string.special", { link = "Special" })
hi(0, "@number", { fg = "#7B9FA2" })
hi(0, "@boolean", { link = "Boolean" })
hi(0, "@constant", { link = "Constant" })
hi(0, "@constant.builtin", { fg = "#7B9FA2", italic = true })
hi(0, "@function", { fg = "#5a93aa" })
hi(0, "@function.call", { fg = "#5a93aa" })
hi(0, "@function.method", { link = "@function.call" })
hi(0, "@function.method.call", { link = "@function.call" })
hi(0, "@function.macro", { link = "@function.call" })
hi(0, "@function.builtin", { fg = "#5a93aa", italic = true })
hi(0, "@keyword", { link = "Keyword" })
hi(0, "@keyword.function", { fg = "#8577af", italic = true })
hi(0, "@keyword.return", { bold = true, fg = "#8577af", italic = true })
hi(0, "@keyword.operator", { bold = true, fg = "#5a93aa", italic = true })
hi(0, "@keyword.conditional.ternary", { bold = true, fg = "#5a93aa" })
hi(0, "@conditional", { link = "Conditional" })
hi(0, "@repeat", { link = "Repeat" })
hi(0, "@exception", { link = "Exception" })
hi(0, "@label", { link = "Label" })
hi(0, "@operator", { link = "Operator" })
hi(0, "@type", { link = "Type" })
hi(0, "@type.builtin", { fg = "#a1cdd8", italic = true })
hi(0, "@module", { fg = "#606267" })
hi(0, "@module.builtin", { fg = "#606267", italic = true })
hi(0, "@variable", { fg = "#e6eaea" })
hi(0, "@variable.builtin", { fg = "#e6eaea", italic = true })
hi(0, "@variable.parameter", { fg = "#B0B0B1" })
hi(0, "@variable.parameter.builtin", { fg = "#B0B0B1", italic = true })
hi(0, "@attribute", { fg = "#afd4de" })
hi(0, "@variable.attribute", { link = "@attribute" })
hi(0, "@attribute.builtin", { fg = "#afd4de", italic = true })
hi(0, "@variable.attribute.builtin", { link = "@attribute.builtin" })
hi(0, "@property", { fg = "#e6eaea" })
hi(0, "@punctuation.bracket", { fg = "#8577af" })
hi(0, "@punctuation.delimiter", { fg = "#8577af" })
hi(0, "@punctuation.special", { fg = "#686377" })
hi(0, "@markup.raw", { fg = "#83AEB0" })
hi(0, "@markup.raw.block", { link = "@markup.raw" })
hi(0, "@markup.link", { fg = "#D6E5E8", underline = true })
hi(0, "@markup.link.url", { fg = "#D6E5E8", underline = true })
hi(0, "@string.special.url", { link = "@markup.link.url" })
hi(0, "@markup.italic", { link = "Italic" })
hi(0, "@markup.strong", { link = "Bold" })
hi(0, "@markup.underline", { link = "Underlined" })
hi(0, "@text.underline", { link = "Underlined" })
hi(0, "@text.strong", { link = "Bold" })
hi(0, "@text.title", { link = "Title" })
hi(0, "@tag", { underline = true })
hi(0, "@tag.builtin", { italic = true, underline = true })
hi(0, "@tag.attribute", { fg = "#afd4de" })
hi(0, "@tag.delimiter", { link = "@comment" })
hi(0, "@tag.html", { underline = false })
hi(0, "@tag.xml", { underline = false })
hi(0, "@conceal", { link = "Comment" })
hi(0, "@conceal.json", { fg = "#7aa4a1" })
hi(0, "@none.python", { link = "Macro" })
hi(0, "@constructor.lua", { link = "@punctuation.bracket" })
hi(0, "@constant.numeric", { link = "Number" })
hi(0, "@include", { link = "PreProc" })

-- Variable is neutral (fg), except in some markup-like languages.
hi(0, "@variable.member", {})
hi(0, "@variable.typst", { fg = "#afd4de" })
hi(0, "@variable.markdown", { fg = "#afd4de" })
hi(0, "@variable.asciidoc", { fg = "#afd4de" })
hi(0, "@variable.latex", { fg = "#afd4de" })
hi(0, "@variable.wgsl", {})

-- Language-specific clears and overrides.
hi(0, "@constant.lua", {})
hi(0, "@constant.python", {})
hi(0, "@type.python", {})
hi(0, "@type.sql", {})
hi(0, "@type.builtin.python", {})
hi(0, "@string.special.symbol.julia", {})
hi(0, "@lsp.type.macro.rust", { link = "@function.call" })
hi(0, "@function.typst", { link = "@function.call" })
hi(0, "@lsp.type.pol.typst", { link = "@variable.typst" })

-- LSP semantic tokens.
hi(0, "@lsp.type.function", { link = "@function.call" })
hi(0, "@lsp.type.method", { link = "@function.call" })
hi(0, "@lsp.type.keyword", { link = "@keyword" })
hi(0, "@lsp.type.operator", { link = "@operator" })
hi(0, "@lsp.type.punct", { link = "Delimiter" })
hi(0, "@lsp.type.string", { link = "@string" })
hi(0, "@lsp.type.escape", { link = "@string.escape" })
hi(0, "@lsp.type.decorator", { link = "PreProc" })
hi(0, "@lsp.type.selfParameter", { link = "@module.builtin" })
hi(0, "@lsp.type.procMacro", { link = "@function.macro" })
hi(0, "@lsp.mod.strong", { bold = true })
hi(0, "@lsp.mod.emph", { italic = true })
hi(0, "@lsp.type.ref", { underline = true })
hi(0, "@lsp.type.link", { underline = true })
hi(0, "@lsp.type.variable.lua", {})

-- Bold for definitions.
hi(0, "@lsp.typemod.function.declaration", { bold = true, fg = "#5a93aa" })
hi(0, "@lsp.typemod.class.definition", { bold = true })
hi(0, "@lsp.typemod.method.definition", { bold = true })
hi(0, "@lsp.typemod.variable.global", { link = "@module.builtin" })
hi(0, "@lsp.typemod.variable.global.lua", { link = "@module.builtin" })
hi(0, "@lsp.typemod.variable.defaultLibrary.lua", { link = "Constant" })
hi(0, "@lsp.type.selfParameter", { link = "@module.builtin" })

-- Completion.
hi(0, "CmpItemAbbrMatch", { bold = true })
hi(0, "CmpItemAbbrMatchFuzzy", { bold = true })

-- Prompt.
hi(0, "Prompt", { fg = "#FF6AC1" })
hi(0, "vimPrompt", { link = "Prompt" })
hi(0, "@punctuation.delimiter.vim", { link = "Prompt" })
hi(0, "NvimPrompt", { link = "Prompt" })
hi(0, "MiniPickPrompt", { link = "Prompt" })
hi(0, "SnacksPickerPrompt", { link = "Prompt" })
hi(0, "FzfLuaFzfPrompt", { link = "Prompt" })
hi(0, "TelescopeSelectionCaret", { link = "Prompt" })

-- Misc.
hi(0, "Quote", { link = "String" })
hi(0, "helpCommand", { link = "@string.special" })
hi(0, "helpExample", { link = "@string.documentation" })
hi(0, "LspSignatureActiveParameter", { link = "Underlined" })
hi(0, "FidgetTitle", { link = "Normal" })
hi(0, "DashboardHeader", { link = "Bold" })
hi(0, "@property.builtin", { underline = true })
hi(0, "@string.builtin", { underline = true })
hi(0, "@string.special.path", { underline = true })
hi(0, "FunctionPath", { fg = "#5a93aa", underline = true })
hi(0, "zshShortDerefPath", { fg = "#5a93aa", underline = true })
hi(0, "Wildcard", { bold = true, fg = "#5a93aa", underline = true })

-- Conflict markers.
hi(0, "ConflictMarkerBegin", { bg = "#1E3236" })
hi(0, "ConflictMarkerOurs", { link = "ConflictMarkerBegin" })
hi(0, "ConflictMarkerCommonAncestors", { bg = "#243A3E" })
hi(0, "ConflictMarkerCommonAncestorsHunk", { link = "ConflictMarkerCommonAncestors" })
hi(0, "ConflictMarkerSeparator", { bg = "#152528", fg = "#e6eaea" })
hi(0, "ConflictMarkerEnd", { bg = "#1B2D32" })
hi(0, "ConflictMarkerTheirs", { link = "ConflictMarkerEnd" })

-- Terminal colours.
vim.g.terminal_color_0 = "#2f3239"
vim.g.terminal_color_1 = "#e85c51"
vim.g.terminal_color_2 = "#7aa4a1"
vim.g.terminal_color_3 = "#fda47f"
vim.g.terminal_color_4 = "#5a93aa"
vim.g.terminal_color_5 = "#8577af"
vim.g.terminal_color_6 = "#a1cdd8"
vim.g.terminal_color_7 = "#ebebeb"
vim.g.terminal_color_8 = "#4e5157"
vim.g.terminal_color_9 = "#eb746b"
vim.g.terminal_color_10 = "#8eb2af"
vim.g.terminal_color_11 = "#fdb292"
vim.g.terminal_color_12 = "#73a3b7"
vim.g.terminal_color_13 = "#FF6AC1"
vim.g.terminal_color_14 = "#afd4de"
vim.g.terminal_color_15 = "#eeeeee"

