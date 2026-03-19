; extends

; Inject miller DSL into single-quoted strings after put/filter/tee verbs
; in mlr commands. Only raw_string (single-quoted) — double-quoted strings
; have zsh variable expansion which conflicts with miller's $ field refs.
(command
  name: (command_name) @_cmd
  argument: (word) @_verb
  .
  argument: (raw_string) @injection.content
  (#eq? @_cmd "mlr")
  (#any-of? @_verb "filter" "put" "tee")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "miller")
  (#set! injection.include-children))
