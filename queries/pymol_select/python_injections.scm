; Inject pymol_select into Python strings in pymol-specific contexts.
; Matches method calls on simple identifiers (cmd.*, pml.*, etc.)
; but NOT chained attributes (sys.path.*, os.path.*).
;
; We capture the full (string) node (not string_content) with
; include-children so that quote delimiters and f-string interpolation
; syntax are included. Combined with injection.combined, the quotes
; act as natural token separators between different strings, preventing
; adjacent identifiers from merging (e.g. "pocket" + "polymer" staying
; separate instead of becoming "pocketpolymer"). The _fallback token
; in the grammar silently consumes the non-pymol characters (", f, {, }).

; obj.method("selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (string) @injection.content
    (#set! injection.language "pymol_select")
    (#set! injection.combined)
    (#set! injection.include-children)))

; obj.method(kw="selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (keyword_argument
      value: (string) @injection.content)
    (#set! injection.language "pymol_select")
    (#set! injection.combined)
    (#set! injection.include-children)))

; obj.method(var + "selection") — string in binary_operator (explicit + concat)
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (binary_operator
      (string) @injection.content)
    (#set! injection.language "pymol_select")
    (#set! injection.combined)
    (#set! injection.include-children)))

; nested: obj.method(a + b + "selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (binary_operator
      (binary_operator
        (string) @injection.content))
    (#set! injection.language "pymol_select")
    (#set! injection.combined)
    (#set! injection.include-children)))

; obj.method(kw=var + "selection")
(call
  function: (attribute
    object: (identifier))
  arguments: (argument_list
    (keyword_argument
      value: (binary_operator
        (string) @injection.content))
    (#set! injection.language "pymol_select")
    (#set! injection.combined)
    (#set! injection.include-children)))

; pymol.cmd.method("selection") — fully qualified
(call
  function: (attribute
    object: (attribute
      object: (identifier) @_pkg
      attribute: (identifier) @_obj))
  arguments: (argument_list
    (string) @injection.content)
  (#eq? @_pkg "pymol")
  (#eq? @_obj "cmd")
  (#set! injection.language "pymol_select")
  (#set! injection.combined)
  (#set! injection.include-children))

; Assignment right-hand side
(assignment
  right: (string) @injection.content
  (#set! injection.language "pymol_select")
  (#set! injection.combined)
  (#set! injection.include-children))

; Concatenated strings in assignments (implicit "a" "b" concat)
(assignment
  right: (concatenated_string
    (string) @injection.content
    (#set! injection.language "pymol_select")
    (#set! injection.combined)
    (#set! injection.include-children)))

; Assignment with string concatenation (explicit + concat)
(assignment
  right: (binary_operator
    (string) @injection.content)
  (#set! injection.language "pymol_select")
  (#set! injection.combined)
  (#set! injection.include-children))
