;extends


"in" @repeat

; fix syntax thinking `x |> f` is piping into a variable.
; NOTE: has to be "lua-match" and not "match", since the latter will use | as OR operator.
(binary_expression
  (_)
  (operator) @_pipe
  (identifier) @function.call
  (#any-of? @_pipe "|>" ".|>"))
; TODO: maybe figure out how to set this as lower priority
; (tried adding (#set! "priority" 0))
; Since piping into types (e.g. x |> String) has the right highlights
; (same as String(x)) but in the wrong order, so the color from builtin type is 
; overruled by function.call.

; support @chain macro recognizing @function.call as well.
(macrocall_expression
  (macro_identifier
    (identifier) @_chain
    (#eq? @_chain "chain"))
  (macro_argument_list
    (identifier) @function.call))

; adding : for `using BSON: @load`
(selected_import
  ":" @punctuation.delimiter)

; highlight the : before a symbol
(quote_expression
  ":" @symbol-prefix)

; (catch_clause
;   ";" @punctuation.delimiter)

