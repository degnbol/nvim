; Inject pymol_select into Python strings in pymol-specific contexts.
; Matches method calls on simple identifiers (cmd.*, pml.*, etc.)
; but NOT chained attributes (sys.path.*, os.path.*).

; obj.method("selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (string (string_content) @injection.content)
    (#set! injection.language "pymol_select")))

; obj.method(kw="selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (keyword_argument
      value: (string (string_content) @injection.content))
    (#set! injection.language "pymol_select")))

; obj.method(var + "selection") — string in binary_operator (explicit + concat)
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (binary_operator
      (string (string_content) @injection.content))
    (#set! injection.language "pymol_select")))

; nested: obj.method(a + b + "selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (binary_operator
      (binary_operator
        (string (string_content) @injection.content)))
    (#set! injection.language "pymol_select")))

; obj.method(kw=var + "selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (keyword_argument
      value: (binary_operator
        (string (string_content) @injection.content)))
    (#set! injection.language "pymol_select")))

; pymol.cmd.method("selection") — fully qualified
(call
  function: (attribute
    object: (attribute
      object: (identifier) @_pkg
      attribute: (identifier) @_obj))
  arguments: (argument_list
    (string (string_content) @injection.content))
  (#eq? @_pkg "pymol")
  (#eq? @_obj "cmd")
  (#set! injection.language "pymol_select"))

; Assignment right-hand side
(assignment
  right: (string (string_content) @injection.content)
  (#set! injection.language "pymol_select"))

; Concatenated strings in assignments (implicit "a" "b" concat)
(assignment
  right: (concatenated_string
    (string (string_content) @injection.content)
    (#set! injection.language "pymol_select")))

; Assignment with string concatenation (explicit + concat)
(assignment
  right: (binary_operator
    (string (string_content) @injection.content))
  (#set! injection.language "pymol_select"))
