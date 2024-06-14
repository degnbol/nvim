;extends

; mark it as error if choice node table contains string instead of a node, e.g. t""
(function_call
 name: (identifier) @_choice
 arguments: (arguments
  (table_constructor
    (field
      value: (string)) @error))
 (#eq? @_choice "c"))

; Mark it as error if t"" contains newline. Use t{"a", "b", "c"} instead.
(function_call
 name: (identifier) @_t
 arguments: (arguments
  (string content: (string_content) @error))
 (#eq? @_t "t")
 (#match? @error "\n"))


