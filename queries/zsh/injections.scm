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

; Inject python into strings after `python3 -c` / `python -c`
; Works with both single-quoted (raw_string) and double-quoted (string) args.
; Two variants: python as the command name, and python as an argument
; (covers `uv run --with pkg python3 -c`, `conda run -n env python3 -c`, etc.)
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "python3" "python")
  (#eq? @_flag "-c")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "python")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-of? @_cmd "python3" "python")
  (#eq? @_flag "-c")
  (#set! injection.language "python")
  (#set! injection.include-children))

(command
  argument: (word) @_py
  .
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_py "python3" "python")
  (#eq? @_flag "-c")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "python")
  (#set! injection.include-children))

(command
  argument: (word) @_py
  .
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-of? @_py "python3" "python")
  (#eq? @_flag "-c")
  (#set! injection.language "python")
  (#set! injection.include-children))

; Inject zsh into `zsh -c '...'` / `zsh -c "..."` / `bash -c ...` / `sh -c ...`
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "zsh" "bash" "sh")
  (#eq? @_flag "-c")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "zsh")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-of? @_cmd "zsh" "bash" "sh")
  (#eq? @_flag "-c")
  (#set! injection.language "zsh")
  (#set! injection.include-children))

; Concatenation form: `zsh -c 'prefix'$var'suffix'`. Each raw_string fragment
; is injected independently; variable_ref siblings are highlighted by the
; outer zsh parser. Fragments may not parse as complete zsh, but most
; highlighting still comes through via error recovery.
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (concatenation
    (raw_string) @injection.content)
  (#any-of? @_cmd "zsh" "bash" "sh")
  (#eq? @_flag "-c")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "zsh")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (concatenation
    (string (string_content) @injection.content))
  (#any-of? @_cmd "zsh" "bash" "sh")
  (#eq? @_flag "-c")
  (#set! injection.language "zsh")
  (#set! injection.include-children))

; Inject julia into `julia -e '...'`
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#eq? @_cmd "julia")
  (#eq? @_flag "-e")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "julia")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#eq? @_cmd "julia")
  (#eq? @_flag "-e")
  (#set! injection.language "julia")
  (#set! injection.include-children))

; Inject R into `Rscript -e '...'` / `R -e '...'`
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "Rscript" "R")
  (#eq? @_flag "-e")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "r")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-of? @_cmd "Rscript" "R")
  (#eq? @_flag "-e")
  (#set! injection.language "r")
  (#set! injection.include-children))

; Inject javascript into `node -e '...'`
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#eq? @_cmd "node")
  (#eq? @_flag "-e")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "javascript")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#eq? @_cmd "node")
  (#eq? @_flag "-e")
  (#set! injection.language "javascript")
  (#set! injection.include-children))

; Inject lua into `lua -e '...'`
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#eq? @_cmd "lua")
  (#eq? @_flag "-e")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#eq? @_cmd "lua")
  (#eq? @_flag "-e")
  (#set! injection.language "lua")
  (#set! injection.include-children))

; Inject awk into the first raw_string argument of awk/gawk/mawk
(command
  name: (command_name) @_cmd
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "awk" "gawk" "mawk")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "awk")
  (#set! injection.include-children))

; Inject lua into `nvim -c "lua ..."` / `nvim -c 'lua ...'`
; The -c flag can appear anywhere in the arg list (after --headless, etc.),
; so no `.` anchor between command name and -c. The `lua ` prefix is skipped
; via #offset! so only the actual Lua code is injected.
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_flag "-c")
  (#lua-match? @injection.content "^'lua[%s]")
  (#offset! @injection.content 0 5 0 -1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_flag "-c")
  (#lua-match? @injection.content "^lua[%s]")
  (#offset! @injection.content 0 4 0 0)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; Inject lua into `nvim +"lua ..."` / `nvim +'lua ...'`
; The + form parses as concatenation(word("+"), string/raw_string).
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (raw_string) @injection.content)
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^'lua[%s]")
  (#offset! @injection.content 0 5 0 -1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (string (string_content) @injection.content))
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^lua[%s]")
  (#offset! @injection.content 0 4 0 0)
  (#set! injection.language "lua")
  (#set! injection.include-children))
