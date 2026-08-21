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
(number) @number

(float) @number.float

(boolean) @boolean
