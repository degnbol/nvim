; Deliberately no (text) capture: a text cell's colour is the injected
; language's business (injections.scm), or nothing.
;
; This file replaces nvim-treesitter's tsv highlights, which give every cell
; `@string`, rather than extending them — an `; extends` file cannot cancel a
; capture. The first non-extending file in runtimepath order is the whole base
; query and every later one is dropped, and $XDG_CONFIG_HOME/nvim precedes the
; site directory that a parser's queries are symlinked into. Not under `after/`
; for the same reason: a non-extending file there sits below the shipped one and
; is ignored.
;
; Side effect once the csv parser is installed: nvim-treesitter's
; queries/csv/highlights.scm is `; inherits: tsv` plus a `","` pattern, and the
; inherited base resolves through the same search — so a csv buffer keeps its
; commas and numbers and loses `@string` on its text cells too, though nothing
; injects there yet (the predicate assumes the tab separator).
;
; The grammar has no comment node. Comments are captured per field and not per
; row, because a row that starts with an empty cell merges into the row above.
; Values are not captured on a comment line, rather than outranked: an
; overlapping highlight keeps every attribute the higher one does not set.
((field) @comment
 (#comment-line? @comment))

((number) @number
 (#not-comment-line? @number))

((float) @number.float
 (#not-comment-line? @number.float))

((boolean) @boolean
 (#not-comment-line? @boolean))
