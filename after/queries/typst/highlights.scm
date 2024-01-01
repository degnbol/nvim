;extends

; ... is captured by @lsp.type.escape.typst that links to comment.
; I don't want ... colored, so I capture it with a higher priority custom 
; capture that gets linked to Normal.
((shorthand) @punct.ellipsis
  (#eq? @punct.ellipsis "...")
  (#set! "priority" 200)); higher priority than 125

