local hl = require "utils/highlights"

hl.def("typstMarkupHeading", "@text.title")
-- Semantic tokens.
-- captures \ at the end of lines, ~ between words and \* that let's you write a literal *. The last one is not ideal but worth it for the other two.
hl.def("@lsp.type.escape", "@comment")
-- undo coloring from above on "..." with help from custom capture in after/queries/typst/highlights.scm
hl.def("@punct.ellipsis", "Normal")
hl.def("@lsp.typemod.punct.emph", "@comment")
hl.def("@lsp.typemod.punct.strong", "@comment")
hl.def("@lsp.type.punct", "Delimiter")
hl.def("@lsp.type.pol", "@variable")
hl.def("@lsp.type.string", "@string")
hl.def("@lsp.type.heading", "@text.title")
hl.def("@lsp.type.number", "@number")

