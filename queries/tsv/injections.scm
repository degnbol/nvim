; Parse the cells of a chemical column as a structure. SMILES and SMARTS are
; both read by the one smarts parser, so the language is constant and the column
; is a predicate — `#chem-column?`, defined in lua/chem/tsv.lua — rather than a
; directive resolving a variable language. A failing predicate drops the whole
; match, so no region is created.
;
; The capture is the `(text)` child and not `(field)`: get_node_ranges masks a
; capture's named children out unless `injection.include-children` is set, and
; `field`'s single child has the same extent, so `(field) @injection.content`
; yields no ranges at all and the region is dropped without a diagnostic.
; Capturing `text` also leaves the numeric cells alone for free.
;
; A cell may be quoted, which is `#unquote!`'s job and not the predicate's.
((field (text) @injection.content)
 (#chem-column? @injection.content)
 (#unquote! @injection.content)
 (#set! injection.language "smarts"))
