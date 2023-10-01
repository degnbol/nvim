;extends
; I want to use bash by default in code blocks in markdown but use the 
; "info_string" if given to set the language.
; So, check fenced_code_block not having an "info_string" subnode.
; Sounds like "negated fields" but I might be confused about field vs child node.
; https://tree-sitter.github.io/tree-sitter/using-parsers#pattern-matching-with-queries
; The solution here works by requiring the specified exact order of child nodes (by using dot), 
; which means there is no other child node allowed between, and that is where 
; the info_string would otherwise go.
(fenced_code_block
  .
  (fenced_code_block_delimiter)
  .
  (block_continuation)
  .
  (code_fence_content) @bash
  )
