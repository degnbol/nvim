((selector) @variable.builtin (#set! priority 200))
((builtin) @variable.builtin (#set! priority 200))
((representation) @type.builtin (#set! priority 200))

((logical_operator) @keyword.operator (#set! priority 200))
((not_operator) @keyword.operator (#set! priority 200))
((expansion_keyword) @keyword.operator (#set! priority 200))
((proximity_keyword) @keyword.operator (#set! priority 200))
((of) @keyword (#set! priority 200))

((comparison_operator) @operator (#set! priority 200))
((wildcard) @operator (#set! priority 200))
((object_ref (object_prefix) @operator) (#set! priority 200))

((number) @number (#set! priority 200))

; macro is a single token — highlight the whole thing as a string escape
; (visually distinct path-like construct within the string)
((macro) @string.escape (#set! priority 200))
(("(" @punctuation.bracket) (#set! priority 200))
((")" @punctuation.bracket) (#set! priority 200))
