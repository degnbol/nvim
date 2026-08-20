; Which element a node is, is its text — there is no per-element node type, so
; per-element colouring belongs in a query layer that matches text. Aromaticity
; is spelled in the letter case, which needs no help from a highlight, so all
; four spellings of an identity share one capture.
[
  (organic_symbol)
  (aromatic_symbol)
  (element_symbol)
  (atomic_number)
] @constant

(wildcard) @string.regexp

; A property of one atom being constrained: stereochemistry, attached hydrogens
; (`H2`), and the predicate primitives (`D3`, `R1`, `v3`, `^3`). Each of the
; latter is a fixed letter naming the property, and its digits are the value —
; a separate node, so the two need no overriding.
(chirality) @attribute
(property_name) @attribute

; Except that the same `H` is the hydrogen atom in the bare spelling, so it takes
; the capture every other identity has, overriding the property above by coming
; after it. RDKit round-trips `[H]` to `[#1]` and `[H+]` to `[#1&+]`, but reads
; `[CH]` as `[C&H1]`, `[HC]` as `[H1&C]` — order is irrelevant — and even
; `[H;+]` as `[H1&+]`, so an explicit conjunction is enough to make it a count.
; The atom is therefore exactly: an optional isotope, `H`, an optional charge, an
; optional map, and nothing else. That is a condition on the bracket's text, not
; on the tree, and matching the bracket says it once. The four spellings differ
; only in how deep (and_high) nests the primitives that follow the `H`.
((bracket_atom (hydrogen (property_name) @constant)) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$"))
((bracket_atom (and_high (hydrogen (property_name) @constant))) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$"))
((bracket_atom (and_high (and_high (hydrogen (property_name) @constant)))) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$"))
((bracket_atom (and_high (and_high (and_high (hydrogen (property_name) @constant))))) @_atom
 (#match? @_atom "^\\[[0-9]*H([-+][0-9]*)*(:[0-9]+)?\\]$"))

; `[#6]` spells an identity in digits, so it sits with the symbols above; what is
; counted goes here. The charge is a signed count, not an arithmetic operator,
; and (count) is every `R1`, `v4`, `H2` digit — overriding the letter's capture
; by coming after it. `^3` is the one that names a hybridisation instead.
(isotope) @number
(charge) @number
(count) @number

; A ring number names a closure bond rather than counting anything, and it is a
; label proper: the digit is free for reuse the moment its pair closes, so
; `C1CC1C1CC1` is two rings, and `%(100)` would name either of them just as well.
; The marker announcing the digit count is spelling, and dims like `[` does.
(ring_number) @label.ring
(ring_bond [
  "%"
  "("
  ")"
] @punctuation.special)

; A map number is the one numeral here that promises stability: within a record
; it names the same atom on both sides of the arrow, which is a reference rather
; than a slot. `:` prefixes it rather than being part of it.
(atom_class ":" @punctuation.delimiter)
(class_number) @variable.atom

(comment) @comment
(name) @string
(header_field) @markup.heading

; Seven of the primitives state an order, which is an electron count — the same
; thing a formal charge states, and a value rather than punctuation. A charge is
; spelled as a number even so, `+` prefixing one (`[C+2]`), where `=` is a
; placeholder standing for the bond as much as for its count: a constant, not a
; numeral. Element identities take that capture too, and the glyph is what tells
; the two apart. Each pattern is scoped to (bond_prim), or it would also claim the
; `:` of an atom class and the `@` of chirality. A dative bond counts the two
; electrons a single bond does and adds which atom both came from, which the arrow
; says by itself.
(bond_prim [
  "-"
  "="
  "#"
  "$"
  ":"
  "->"
  "<-"
] @constant)

; Any bond, as `*` is any atom.
(bond_prim "~" @string.regexp)

; The rest constrain a bond without counting anything, and take the capture their
; atom-side counterpart has: `@` is ring membership, as the `R` primitive is for
; an atom, and `/` and `\` are cis/trans direction, as chirality is.
(bond_prim [
  "@"
  "/"
  "\\"
  "\\\\"
] @attribute)

; The exception: `.` separates two components, so it is the one bond glyph that
; counts no electrons. Scoped, as `.` also delimits the extended layer's fields.
(chain "." @punctuation.special)

; Negation and the two conjunctions, `&` binding tighter than `;`. All symbolic,
; and `!` prefixes bond primitives too. `;` is scoped for the same reason as `.`.
[
  "!"
  "&"
] @operator

(and_low ";" @operator)
(bond ";" @operator)

; Disjunction reads as a listing of what is allowed at one position — `[C,N]` —
; so it delimits rather than operates.
"," @punctuation.delimiter

; Grouping, a branch of the structure or a nested pattern. Scoped, or these would
; also claim the parens of a long ring number, which are spelling.
(branch ["(" ")"] @punctuation.bracket)
(recursive ["(" ")"] @punctuation.bracket)

; Around annotations of a single atom — `[13CH3+]` — the brackets are a spelling
; requirement carrying nothing themselves, so they dim like the `.` above.
(bracket_atom ["[" "]"] @punctuation.special)

; They delimit something once the atom is a listing of alternatives, and then
; read as brackets. Only `,` makes a listing: `&` and `;` conjoin constraints on
; the one atom and leave it a single specification. Matching the text rather than
; the child node type finds a `,` at any nesting depth, and is exact: no
; primitive token contains one. Order does the overriding, later pattern winning.
((bracket_atom ["[" "]"] @punctuation.bracket) @_list
 (#match? @_list ","))

; The recursion marker, turning a nested pattern into one atom primitive. Shares
; its glyph with the quadruple bond above, so the two must not share a colour.
(recursive "$" @operator)

; The one reserved structural token of a record: it separates a reaction into
; reactant, agent and product sides. Emphasised further in this config's
; after/ layer, where reading left from right is worth a loud divider. Scoped,
; because a `>` also occurs as the query operator of a data sgroup below, and
; the dative bonds `->` and `<-` are single tokens of their own.
(structure ">" @punctuation.delimiter)

; The extended layer: what annotates a structure rather than states it. Its
; punctuation reuses glyphs the structure spells other things with — `;` conjoins
; inside a bracket atom, `.` separates two components — so each pattern here is
; scoped to the field whose delimiter it is.
(cx_layer "|" @punctuation.bracket)

; Real coordinates: `,` between the components of one point, `;` between points.
(cx_coordinates ["(" ")"] @punctuation.bracket)
(cx_coordinates ";" @punctuation.delimiter)
(coordinate) @number.float

; Per-atom text, one entry per atom in writing order: a label naming the atom, or
; in the `$_AV:` spelling the molfile value of it. Both are free text.
(cx_atom_labels "$" @punctuation.bracket)
(cx_atom_labels ";" @punctuation.delimiter)
(cx_atom_values "$_AV:" @keyword)
(cx_atom_values "$" @punctuation.bracket)
(cx_atom_values ";" @punctuation.delimiter)
(label) @string
(escape) @string.escape

; `atomProp:0.name.value` is the one field carrying a key of the writer's own
; choosing, so the key is a property and the value the string it is set to. The
; marker is not that key: it announces that keys follow.
(cx_atom_props "atomProp" @keyword)
(cx_atom_props ":" @punctuation.delimiter)
(atom_prop "." @punctuation.delimiter)
(prop_name) @property
(prop_value) @string

; Substance groups. Their text fields are positional — name, data, query
; operator, unit, tag — so none of them names a key, while the type code comes
; from a closed vocabulary and names a kind of group.
(cx_data_sgroup "SgD" @keyword)
(cx_data_sgroup ":" @punctuation.delimiter)
(cx_polymer_sgroup "Sg" @keyword)
(cx_polymer_sgroup ":" @punctuation.delimiter)
(sgroup_type) @type
(sgroup_field) @string

; Every field marker comes from the same closed vocabulary the parser dispatches
; on, and each states which kind of annotation the indices after it are — so they
; are the reserved words of this layer, whatever shape their arguments take.
; Unrecognised ones are skipped by RDKit rather than read as a name of any kind.
(cx_name) @keyword
(cx_field [":" "."] @punctuation.delimiter)
(index) @number
