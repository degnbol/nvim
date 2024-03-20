;extends

"in" @repeat

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

; more focus to the dot that indicates a broadcast
(broadcast_call_expression
  "." @punctuation.special)

