;extends
; temp solution, highlight inline code with single color
(code_span) @code_span

; Un-conceal tildes inside strikethrough nodes by default — the upstream
; parser pairs unrelated single tildes (e.g. `~14 vs ~7`) as strikethrough
; delimiters, and we want those visible as plain text.
; https://github.com/tree-sitter-grammars/tree-sitter-markdown/issues/236
((strikethrough
  (emphasis_delimiter) @conceal)
  (#set! priority 200)
  (#set! conceal "~"))

; Re-conceal tildes for the genuine `~~double~~` case, which parses as a
; strikethrough containing a nested strikethrough — never produced by stray
; single tildes. Higher priority than the un-conceal above so nested wins.
; Tree-sitter sibling patterns are order-sensitive: leading and trailing
; outer tildes need separate patterns.
(strikethrough
  (emphasis_delimiter) @conceal
  (strikethrough)
  (#set! priority 250)
  (#set! conceal ""))

(strikethrough
  (strikethrough)
  (emphasis_delimiter) @conceal
  (#set! priority 250)
  (#set! conceal ""))

(strikethrough
  (strikethrough
    (emphasis_delimiter) @conceal)
  (#set! priority 250)
  (#set! conceal ""))

; Highlight only the inner of a nested (double-tilde) strikethrough as struck
; text. The base @markup.strikethrough capture from upstream is cleared in
; ftplugin/markdown.lua to suppress single-tilde false positives, so this
; capture (linked to a separate group) is the only way struck text renders.
(strikethrough
  (strikethrough) @markup.strikethrough.double)

; Highlight the angle brackets of autolinks as @punctuation.bracket. The
; grammar parses `<https://...>` and `<email@...>` as a single leaf token,
; so #head!/#tail! narrow the capture range to the first/last byte. Priority
; above the default 100 so the bracket style wins over the wider
; @markup.link.url that covers the whole token.
((uri_autolink) @punctuation.bracket
  (#head! @punctuation.bracket 1)
  (#set! priority 110))
((uri_autolink) @punctuation.bracket
  (#tail! @punctuation.bracket 1)
  (#set! priority 110))
((email_autolink) @punctuation.bracket
  (#head! @punctuation.bracket 1)
  (#set! priority 110))
((email_autolink) @punctuation.bracket
  (#tail! @punctuation.bracket 1)
  (#set! priority 110))

; Same for inline-link punctuation `[label](url)` and image `![alt](url)`.
; Visible at conceallevel=0; the upstream `(#set! conceal "")` still hides
; them at level >= 1 (we don't override the conceal, only add a highlight).
(inline_link
  [
    "["
    "]"
    "("
    ")"
  ] @punctuation.bracket
  (#set! priority 110))

(image
  [
    "!"
    "["
    "]"
    "("
    ")"
  ] @punctuation.bracket
  (#set! priority 110))

