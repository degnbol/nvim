; Every node whose text is an element identity. Not `highlights.scm` and no
; `; extends`, so it stays out of the highlights chain: which element a node is
; cannot be a capture name without a pattern per element, and the highlighter
; reads the group from the capture name alone. lua/chem/highlight.lua looks the
; text up instead and sets the group on an extmark.
;
; On disk rather than built in Lua so `:InspectTree` and `:EditQuery` can be
; pointed at it with no Lua having run.
[
  (organic_symbol)
  (aromatic_symbol)
  (element_symbol)
  (atomic_number)
] @atom
