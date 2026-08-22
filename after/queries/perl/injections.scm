;extends

; -----------------------------------------------------------------------------
; Inject regex into every pattern: `m//`, `s///` and `qr//` all node theirs as
; (regexp_content), so one pattern covers all three. `tr///` uses
; (transliteration_content) — a character map, not a regex — and is untouched.
; The node spans the pattern only, delimiters excluded, so no `#trim!` is needed.
;
; The base query paints the whole pattern flat @string.regexp; the injected layer
; renders over it (deeper layer wins at equal priority), so classes, groups and
; quantifiers get their own colours.
;
; include-children keeps the pattern contiguous: escapes and interpolations are
; named children of (regexp_content), and without it they would be cut out of the
; region, leaving the regex parser `[` + `.-]` for `[\w\/.-]`. Interpolated
; variables keep their perl colours via after/queries/perl/highlights.scm.
;
; The grammar models PCRE/ECMAScript — perl's own dialect, so the BRE/ERE
; divergences `grammar_models_flavour` (lua/utils/treesitter.lua) gates the grep
; patterns on do not arise here. What perl has and the grammar lacks —
; possessive quantifiers (`a++`), recursion (`(?1)`), inline comments (`(?#…)`),
; code blocks (`(?{…})`), branch reset (`(?|…)`) — parses to an ERROR node, so
; those patterns degrade to partial colours over the base @string.regexp rather
; than mis-colouring. `\K`, `\A`, `\z`, `\h`, `\R` render as inert identity
; escapes, which is no worse than the flat colour they get without the injection.
;
; Under the `/x` modifier free whitespace and `#` comments are insignificant to
; perl but literal to the regex parser, so such a comment takes pattern colours.
; Not worth gating on: the content is still a regex, unlike an `/e` replacement
; (which is perl code — the base query injects perl for it).
; -----------------------------------------------------------------------------
((regexp_content) @injection.content
  (#set! injection.language "regex")
  (#set! injection.include-children))
