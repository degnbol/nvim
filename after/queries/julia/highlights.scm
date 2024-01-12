;extends

"in" @repeat

; fix syntax thinking `x |> f` is piping into a variable.
; NOTE: has to be "lua-match" and not "match", since the latter will use | as OR operator.
(binary_expression
  (_)
  (operator) @_pipe
  (identifier) @function.call
  (#lua-match? @_pipe "\.?|>"))

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

