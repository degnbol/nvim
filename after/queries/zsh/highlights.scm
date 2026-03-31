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
