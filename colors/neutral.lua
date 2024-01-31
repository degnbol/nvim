local hl = require "utils/highlights"

-- mods to default nvim colorscheme
vim.cmd "hi clear"

hl.set("Keyword", {fg="#9a808f", italic=true})
hl.link("@conditional", "Keyword")
hl.link("@repeat", "Keyword")
hl.link("@include", "Keyword")
hl.link("Exception", "Keyword")
hl.set("Type", {fg="#9a808f", bold=true})
hl.link("Operator", "Type")
hl.set("@keyword.function", {fg="#9a808f"})
hl.set("Preproc", {fg="#8888ff"})
hl.link("Identifier", "@variable")
hl.link("Special", "PreProc")
hl.link("Number", "String")
hl.link("Constant", "String")
hl.link("@constant.builtin", "String")
