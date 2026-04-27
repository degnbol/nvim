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

; Inject vim (vimscript) into `nvim -c '...'` / `nvim -c "..."` and
; `nvim +'...'` / `nvim +"..."`. The vim parser's own injections.scm handles
; nested lua/python/ruby for `:lua print(1)`, `:python << EOF ... EOF`, etc.
;
; Special case: multi-line `nvim -c 'lua\nCODE\n'` is NOT valid vim syntax
; (vim's heredoc form is `:lua << EOF\n...\nEOF`), but it's a common shell
; shorthand. A separate pattern below injects lua directly for that form.

; -c '...'  → vim
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_flag "-c")
  (#not-lua-match? @injection.content "^'lua\n")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "vim"))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_flag "-c")
  (#not-lua-match? @injection.content "^lua\n")
  (#set! injection.language "vim"))

; +'...' / +"..." form: concatenation(word("+"), string/raw_string).
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (raw_string) @injection.content)
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#not-lua-match? @injection.content "^'lua\n")
  (#offset! @injection.content 0 1 0 -1)
  (#set! injection.language "vim"))

(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (string (string_content) @injection.content))
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#not-lua-match? @injection.content "^lua\n")
  (#set! injection.language "vim"))

; Multi-line lua: `nvim -c 'lua\nCODE\n'` and `nvim +"lua\nCODE\n"` etc.
; #trim! is a custom directive (see plugin/treesitter.lua) that skips a
; byte-prefix and byte-suffix from the captured node, computing a (row, col,
; byte) range with all three coordinates consistent. `#offset!` can't do this
; because it does naive (row+drow, col+dcol) arithmetic and the col delta
; needed to reach column 0 of the next line varies with the surrounding text.
;
; raw_string `'lua\n...\n'` — skip 5 bytes (`'lua\n`), strip 1 byte (`'`).
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_flag "-c")
  (#lua-match? @injection.content "^'lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua"))

; double-quoted `"lua\n...\n"` — tree-sitter-zsh splits this into multiple
; `string_content` nodes with `"`s as siblings, so capturing string_content
; directly only sees one line at a time. Match the whole `string` node and
; trim 5 bytes (`"lua\n`) from start, 1 byte (`"`) from end.
; `include-children` keeps the metadata.range from #trim! intact (otherwise
; the range gets masked by the inner string_content/`"` children).
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string) @injection.content
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_flag "-c")
  (#lua-match? @injection.content "^\"lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; +'lua\n...\n' concatenation form
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (raw_string) @injection.content)
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^'lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua"))

(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (string) @injection.content)
  (#any-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^\"lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; -----------------------------------------------------------------------------
; Heredoc body injection based on file-redirect extension
;
; `cat > foo.lua <<EOF ... EOF` — file_redirect sibling of heredoc_redirect.
; `cat <<EOF > foo.lua ... EOF` — file_redirect nested inside heredoc_redirect.
; Destination may be a bare word, double-quoted string, or raw_string — match
; any node type via `(_)` and let the regex handle trailing quotes.
; The base zsh query already uses `heredoc_end` as the language, so
; `<<LUA ... LUA` still works without an extension hint.
; -----------------------------------------------------------------------------

; lua
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.lua[\"']?$")
  (#set! injection.language "lua"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.lua[\"']?$")
  (#set! injection.language "lua"))

; python
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.py[\"']?$")
  (#set! injection.language "python"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.py[\"']?$")
  (#set! injection.language "python"))

; julia
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.jl[\"']?$")
  (#set! injection.language "julia"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.jl[\"']?$")
  (#set! injection.language "julia"))

; R
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.[rR][\"']?$")
  (#set! injection.language "r"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.[rR][\"']?$")
  (#set! injection.language "r"))

; javascript (.js, .mjs, .cjs)
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.[cm]?js[\"']?$")
  (#set! injection.language "javascript"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.[cm]?js[\"']?$")
  (#set! injection.language "javascript"))

; zsh
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.zsh[\"']?$")
  (#set! injection.language "zsh"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.zsh[\"']?$")
  (#set! injection.language "zsh"))

; bash
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.bash[\"']?$")
  (#set! injection.language "bash"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.bash[\"']?$")
  (#set! injection.language "bash"))

; sh (.sh → bash parser; literal `.` prevents matching .zsh / .bash)
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#lua-match? @_dest "%.sh[\"']?$")
  (#set! injection.language "bash"))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#lua-match? @_dest "%.sh[\"']?$")
  (#set! injection.language "bash"))
