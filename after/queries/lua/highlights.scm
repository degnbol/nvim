;extends

; mark it as error if choice node table contains string instead of a node, e.g. t""
(function_call
 name: (identifier) @_choice
 arguments: (arguments
  (table_constructor
    (field
      value: (string)) @error))
 (#eq? @_choice "c"))

