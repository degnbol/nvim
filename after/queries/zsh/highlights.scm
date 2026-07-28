;extends

; --- numbers ---

; The base query only matches integers. Capture decimals like 0.5, 3.14, .5
((word) @number
  (#lua-match? @number "^[0-9]*%.[0-9]+$"))

; --- paths ---

; Command names that are paths: $ROOT/src/script.py, ./run.sh
; Path on the whole concatenation (underline), function.call on the inner word
; to override the base query's (concatenation (word) @string).
(command_name
  (concatenation) @string.special.path
  (#lua-match? @string.special.path "/"))

(command_name
  (concatenation
    (word) @function.call
    (#lua-match? @function.call "/")))

(command_name
  (word) @string.special.path
  (#lua-match? @string.special.path "/"))

(command_name
  (word) @function.call
  (#lua-match? @function.call "/"))

; Words/concatenations containing / are paths (command args, for-loop items, etc.)
(command
  argument: (word) @string.special.path
  (#lua-match? @string.special.path "/"))

(command
  argument: (concatenation) @string.special.path
  (#lua-match? @string.special.path "/"))

; Filenames with extensions: word.ext (letter after dot, avoids 0.5)
(command
  argument: (word) @string.special.path
  (#lua-match? @string.special.path "[a-zA-Z0-9_*-]+%.[a-zA-Z][a-zA-Z0-9._-]*$"))

(command
  argument: (concatenation) @string.special.path
  (#lua-match? @string.special.path "%.[a-zA-Z][a-zA-Z0-9._-]*$"))

; Dotfiles: .gitignore, .bashrc, etc.
(command
  argument: (word) @string.special.path
  (#lua-match? @string.special.path "^%.[a-zA-Z_]"))

; Override base (concatenation (word) @string) when the word has a file extension
; or starts with / (path component). Covers ${fname:r}_pocket.tsv, $ROOT/src/...
(concatenation
  (word) @string.special.path
  (#lua-match? @string.special.path "%.[a-zA-Z][a-zA-Z0-9._-]*$"))

(concatenation
  (word) @string.special.path
  (#lua-match? @string.special.path "/"))

; Paths in variable assignments: ROOT=./path or VAR=/some/path
(variable_assignment
  (word) @string.special.path
  (#lua-match? @string.special.path "/"))

(variable_assignment
  (concatenation) @string.special.path
  (#lua-match? @string.special.path "/"))

; for-loop items: `for f in file.tsv other.csv`
(for_statement
  (word) @string.special.path
  (#lua-match? @string.special.path "/"))

(for_statement
  (word) @string.special.path
  (#lua-match? @string.special.path "[a-zA-Z0-9_*-]+%.[a-zA-Z][a-zA-Z0-9._-]*$"))

; --- $expansions inside injected regions ---

; Keep shell $expansions in their zsh colours when they fall inside an injected
; region — interpreter double-quoted strings and heredoc bodies (see
; queries/zsh/injections.scm). Those injections cover the whole string/body
; including the $expansion bytes, so the injected language would otherwise
; colour them. Re-assert the base zsh captures at priority 101 (> injected
; default 100) so the outer language wins for those cells. This reuses the base
; groups and mirrors the base query's split — `$`/`${`/`}` sigils =
; @punctuation.special, name = @variable — rather than flattening the whole node
; to one @variable. Scoped to string/heredoc_body so ordinary shell code keeps
; the base @variable.builtin/@constant refinements. Command substitutions
; ($(...)) are left to their own highlighting.
((string
   (variable_ref "$" @punctuation.special (simple_variable_name) @variable))
 (#set! priority 101))
((string
   (expansion "${" @punctuation.special
     (simple_variable_name) @variable "}" @punctuation.special))
 (#set! priority 101))
((heredoc_body
   (variable_ref "$" @punctuation.special (simple_variable_name) @variable))
 (#set! priority 101))
((heredoc_body
   (expansion "${" @punctuation.special
     (simple_variable_name) @variable "}" @punctuation.special))
 (#set! priority 101))
