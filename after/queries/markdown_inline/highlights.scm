;extends
; temp solution, highlight inline code with single color
(code_span) @code_span

; Un-conceal tildes inside strikethrough nodes — the upstream parser
; incorrectly pairs single tildes as strikethrough delimiters.
; https://github.com/tree-sitter-grammars/tree-sitter-markdown/issues/236
((strikethrough
  (emphasis_delimiter) @conceal)
  (#set! priority 200)
  (#set! conceal "~"))

