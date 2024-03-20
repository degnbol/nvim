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
  ":" @punctuation)

; more focus to the dot that indicates a broadcast.
; Using @operator as opposed to e.g. @punctuation since operator is used for .|> and for julia regex syntax.
(broadcast_call_expression
  "." @operator)

