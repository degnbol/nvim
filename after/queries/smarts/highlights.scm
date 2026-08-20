; extends

; The style layer. Element identity is node text rather than node type, so it
; cannot be a capture name at all: lua/chem/highlight.lua looks the text up and
; sets the group on an extmark. Every pattern here is one no enumeration of the
; periodic table could produce. Priority over treesitter's default 100 rather
; than relying on pattern order (:help treesitter-highlight-priority).

; Hydrogen is the one element with no symbol node: the same `H` is the atom in the
; bare spelling and a count of attached hydrogens everywhere else, so the two
; readings take different colours. Which is which is a condition on the bracket's
; whole text — see the grammar's own query, where the same predicate says why.
; The `[#1]` spelling is never ambiguous, and falls out of the text lookup.
((bracket_atom (hydrogen (property_name) @chem.element.hydrogen)) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$")
 (#set! priority 101))
((bracket_atom (and_high (hydrogen (property_name) @chem.element.hydrogen))) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$")
 (#set! priority 101))
((bracket_atom (and_high (and_high (hydrogen (property_name) @chem.element.hydrogen)))) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$")
 (#set! priority 101))
((bracket_atom (and_high (and_high (and_high (hydrogen (property_name) @chem.element.hydrogen))))) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$")
 (#set! priority 101))

; Two reaction arrows per line at most, and which side an atom is on decides its
; meaning, so they get a divider rather than a punctuation colour. Scoped to the
; structure, since a `>` also occurs in the extended layer as a query operator.
((structure ">" @chem.reaction)
 (#set! priority 101))

; Weight marks a bond whatever it states, leaving the hue to say which of an
; order, a topology, a direction or a wildcard that is. It also tells the orders
; apart from a charge or an element symbol, which share their colour by design.
; The `!` of `!-` is an operator, not part of the bond.
((bond_prim _ @chem.bond)
 (#not-eq? @chem.bond "!")
 (#set! priority 101))

; The atom a recursive pattern is anchored to: "a SMARTS starting with the atom
; of interest" (Daylight theory manual §4.4). It is the only atom of the nested
; pattern that ends up in the match, and which one it is depends on the writing
; order alone — `[$(*O)]` is the carbon, `[$(O*)]` the oxygen. Nested chains are
; not direct children of this one, so the first atom is the only one matched.
((recursive (chain (branched_atom [
  (organic_symbol)
  (aromatic_symbol)
  (wildcard)
  (bracket_atom)
] @chem.anchor)))
 (#set! priority 101))
