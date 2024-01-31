local hi = require "utils/highlights"

hi.hide("qfSeparatorLeft")
hi.hide("qfSeparatorRight")

hi.def("qfFileName", "Directory")
hi.def("qfLineNr", "LineNr")
hi.def("qfError", "DiagnosticError")
hi.def("qfWarning", "DiagnosticWarn")
hi.def("qfInfo", "DiagnosticInfo")
hi.def("qfNote", "DiagnosticHint")

