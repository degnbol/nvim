(selector) @variable.builtin
(builtin) @variable.builtin
(representation) @type.builtin

(logical_operator) @keyword.operator
(not_operator) @keyword.operator
(expansion_keyword) @keyword.operator
(proximity_keyword) @keyword.operator
(of) @keyword

(comparison_operator) @operator
(wildcard) @operator
(object_ref (object_prefix) @operator)

(number) @number

; macro is a single token — highlight the whole thing as a string escape
; (visually distinct path-like construct within the string)
(macro) @string.escape
"(" @punctuation.bracket
")" @punctuation.bracket
