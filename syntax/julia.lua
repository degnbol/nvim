local hl = require "utils/highlights"

-- Syntax highlights provide JuliaSemicolon, JuliaOperator, etc.
-- We currently enable it since treesitter doesn't capture ; as expression delimiter.
hl.def("JuliaSemicolon", "@punctuation.delimiter")
-- We disable JuliaOperator since it captures ...
-- Syntax highlights are too crude and doesn't have punctuation as a concept.
-- We prefer @punctuation.special for ... since operator will highlight the dots like a broadcast dot, and they are very different.
hl.def("JuliaOperator", "Normal")
