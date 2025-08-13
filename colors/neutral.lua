local hi = require "utils/highlights"

-- mods to default nvim colorscheme
vim.cmd "hi clear"

hi.set("Keyword", { fg = "#9a808f", italic = true })
hi.link("@conditional", "Keyword")
hi.link("@repeat", "Keyword")
hi.link("@include", "Keyword")
hi.link("Exception", "Keyword")
hi.set("Type", { fg = "#9a808f", bold = true })
hi.link("Operator", "Type")
hi.set("@keyword.function", { fg = "#9a808f" })
hi.set("Preproc", { fg = "#8888ff" })
hi.link("Identifier", "@variable")
hi.link("Special", "PreProc")
hi.link("Number", "String")
hi.link("Constant", "String")
hi.link("@constant.builtin", "String")
