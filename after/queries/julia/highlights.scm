;extends
"in" @repeat

; fix syntax thinking `x |> f` is piping into a variable.
(binary_expression
  (_)
  (operator) @pipe
  (identifier) @function.call
  (#lua-match? @pipe "|>"))
; NOTE: has to be lua-match and not match, since the latter will use | as OR operator.

