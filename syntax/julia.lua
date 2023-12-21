local hl = require "utils/highlights"

-- @symbol-prefix defined in after/query/julia/highlights.scm
hl.link("@symbol", "String")
hl.link("@symbol-prefix", "@variable")

