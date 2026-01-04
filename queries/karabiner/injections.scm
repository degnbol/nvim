; Inject bash highlighting into shell_command string values
; Requires bash treesitter parser to be installed
((pair
  key: (string (string_content) @_key)
  value: (string (string_content) @injection.content))
  (#eq? @_key "shell_command")
  (#set! injection.language "bash"))
