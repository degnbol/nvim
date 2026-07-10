; extends

; Inject miller DSL into single-quoted strings after put/filter/tee verbs
; in mlr commands. Only raw_string (single-quoted) — double-quoted strings
; have zsh variable expansion which conflicts with miller's $ field refs.
(command
  name: (command_name) @_cmd
  argument: (word) @_verb
  .
  argument: (raw_string) @injection.content
  (#any-basename-of? @_cmd "mlr")
  (#any-of? @_verb "filter" "put" "tee")
  (#trim! @injection.content 1 1)
  (#set! injection.language "miller")
  (#set! injection.include-children))

; A verb chained via `+\` (backslash immediately after `+`, no space) makes
; tree-sitter-zsh split the chain: `mlr … +` becomes one command and the verb
; starts a new one, so `put`/`filter` land as command_name rather than an
; argument to `mlr`. Inject miller when such a verb command_name is immediately
; followed by a single-quoted string. `tee` is excluded here (unlike above): it
; is a real command whose miller verb takes a filename, not DSL, so a bare
; `tee '…'` command would be a false positive.
(command
  name: (command_name) @_verb
  .
  argument: (raw_string) @injection.content
  (#any-of? @_verb "filter" "put")
  (#trim! @injection.content 1 1)
  (#set! injection.language "miller")
  (#set! injection.include-children))

; -----------------------------------------------------------------------------
; Interpreter `<flag> '<code>'` injections (python -c, zsh -c, julia -e, ...).
;
; #inject-interp! (see lua/utils/treesitter.lua) resolves the interpreter and
; target language in Lua: it gates @_flag to a single-dash cluster ending in a
; command char (`-c`/`-e`, plus clusters like `-lc`), walks back past other
; flags to the interpreter token, and looks it up in a basename → language
; table. This is O(1) in concurrent partial matches regardless of command
; length — a floating `@_interp` capture is O(command length) and silently
; starves the injection out past tree-sitter's match_limit on long commands
; (see notes/PLAN-zsh-injection-match-explosion.md).
;
; There is no static injection.language: an off-table interpreter or a
; non-command flag leaves it unset, so the capture is ignored (same contract as
; #inject-by-ext!). That lets `nvim -c` / `sqlite3 db 'sql'` / `grep -e` fall
; through to their own patterns or to no injection.
;
; raw_string (single-quoted) and string (double-quoted) need separate patterns
; (differing #trim!/string_content handling), plus the `'a'$var'b'`
; concatenation form.
; -----------------------------------------------------------------------------
(command
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#inject-interp! @_flag)
  (#trim! @injection.content 1 1)
  (#set! injection.include-children))

(command
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#inject-interp! @_flag)
  (#set! injection.include-children))

; Concatenation form: `zsh -c 'prefix'$var'suffix'`. Each raw_string / string
; fragment is injected independently; variable_ref siblings are highlighted by
; the outer zsh parser. Fragments may not parse as complete code, but most
; highlighting still comes through via error recovery.
(command
  argument: (word) @_flag
  .
  argument: (concatenation
    (raw_string) @injection.content)
  (#inject-interp! @_flag)
  (#trim! @injection.content 1 1)
  (#set! injection.include-children))

(command
  argument: (word) @_flag
  .
  argument: (concatenation
    (string (string_content) @injection.content))
  (#inject-interp! @_flag)
  (#set! injection.include-children))

; Inject awk into the first raw_string argument of awk/gawk/mawk
(command
  name: (command_name) @_cmd
  argument: (raw_string) @injection.content
  (#any-basename-of? @_cmd "awk" "gawk" "mawk")
  (#trim! @injection.content 1 1)
  (#set! injection.language "awk")
  (#set! injection.include-children))

; Inject jq into every raw_string / string argument of jq / gojq. The filter
; is whichever quoted arg(s) the user passes; arg-bearing flags like `--arg`,
; `--argjson` also take quoted values, but those parse as jq string literals
; without breaking, so highlighting them as jq is harmless.
(command
  name: (command_name) @_cmd
  argument: (raw_string) @injection.content
  (#any-basename-of? @_cmd "jq" "gojq")
  (#trim! @injection.content 1 1)
  (#set! injection.language "jq")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (string (string_content) @injection.content)
  (#any-basename-of? @_cmd "jq" "gojq")
  (#set! injection.language "jq")
  (#set! injection.include-children))

; Inject vim (vimscript) into `nvim -c '...'`, `nvim --cmd '...'`, and
; `nvim +'...'` (and the double-quoted equivalents). The vim parser's own
; injections.scm handles nested lua/python/ruby for `:lua print(1)`,
; `:python << EOF ... EOF`, etc.
;
; Special case: multi-line `nvim -c 'lua\nCODE\n'` is NOT valid vim syntax
; (vim's heredoc form is `:lua << EOF\n...\nEOF`), but it's a common shell
; shorthand. A separate pattern below injects lua directly for that form.

; -c / --cmd '...'  → vim
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-basename-of? @_cmd "nvim" "vim")
  (#any-of? @_flag "-c" "--cmd")
  (#not-lua-match? @injection.content "^'lua\n")
  (#trim! @injection.content 1 1)
  (#set! injection.language "vim"))

(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: (string (string_content) @injection.content)
  (#any-basename-of? @_cmd "nvim" "vim")
  (#any-of? @_flag "-c" "--cmd")
  (#not-lua-match? @injection.content "^lua\n")
  (#set! injection.language "vim"))

; +'...' / +"..." form: concatenation(word("+"), string/raw_string).
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (raw_string) @injection.content)
  (#any-basename-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#not-lua-match? @injection.content "^'lua\n")
  (#trim! @injection.content 1 1)
  (#set! injection.language "vim"))

(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (string (string_content) @injection.content))
  (#any-basename-of? @_cmd "nvim" "vim")
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
  (#any-basename-of? @_cmd "nvim" "vim")
  (#any-of? @_flag "-c" "--cmd")
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
  (#any-basename-of? @_cmd "nvim" "vim")
  (#any-of? @_flag "-c" "--cmd")
  (#lua-match? @injection.content "^\"lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; +'lua\n...\n' concatenation form
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (raw_string) @injection.content)
  (#any-basename-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^'lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua"))

(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus (string) @injection.content)
  (#any-basename-of? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^\"lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; -----------------------------------------------------------------------------
; Heredoc body injection by file-redirect extension. #inject-by-ext! (see
; lua/utils/treesitter.lua) maps the destination's extension → filetype →
; parser language at runtime, so any filetype neovim recognises with an
; installed parser works (lua, python, typst, json, …) — no per-language rule.
;
; `cat > foo.lua <<EOF ... EOF` — file_redirect sibling of heredoc_redirect.
; `cat <<EOF > foo.lua ... EOF` — file_redirect nested inside heredoc_redirect.
; Destination may be a bare word, double-quoted string, or raw_string — match
; any node type via `(_)`; the directive strips trailing quotes. Unknown
; extensions leave the language unset, so the base zsh query's heredoc_end-tag
; injection (`<<LUA ... LUA`) still applies.
; -----------------------------------------------------------------------------
(redirected_statement
  (file_redirect destination: (_) @_dest)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#inject-by-ext! @_dest))

(heredoc_redirect
  (file_redirect destination: (_) @_dest)
  (heredoc_body) @injection.content
  (#inject-by-ext! @_dest))

; -----------------------------------------------------------------------------
; `nvim --headless -l /dev/stdin <<EOF ... EOF` — inject lua into the heredoc
; body. `nvim -l <file>` executes <file> as a lua script; with `/dev/stdin`
; (or `-`) the heredoc body is what runs. Anchor `-l` adjacent to its value
; only — preceding args (`--headless`, `-u NONE`, …) are unconstrained.
; A trailing `2>&1` parses as a file_redirect nested inside heredoc_redirect;
; its destination is a (number), which `vim.filetype.match` can't resolve, so
; the earlier `#inject-by-ext!` rule is a no-op on it.
; -----------------------------------------------------------------------------
(redirected_statement
  body: (command
    name: (command_name) @_cmd
    argument: (word) @_lflag
    .
    argument: (word) @_stdin)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#any-basename-of? @_cmd "nvim" "vim")
  (#eq? @_lflag "-l")
  (#any-of? @_stdin "/dev/stdin" "-")
  (#set! injection.language "lua"))

; -----------------------------------------------------------------------------
; Inject SQL into the query argument of `sqlite3 <db> '<sql>'`.
;
; sqlite3's invocation is `sqlite3 [flags...] <database> <sql>` — the SQL is
; always the LAST argument, preceded by the database path (e.g.
; `"file:zotero.sqlite?immutable=1"`). Matching the last argument (trailing `.`)
; that is itself preceded by another argument injects only the SQL and never the
; db path or a lone `sqlite3 <db>` (interactive) invocation. The preceding
; `argument: (_)` must NOT be anchored adjacent — an anchor binds it to the
; first argument and fails when flags (`-readonly`, `-json`, …) precede the db.
; -----------------------------------------------------------------------------
(command
  name: (command_name) @_cmd
  argument: (_)
  argument: (raw_string) @injection.content
  .
  (#any-basename-of? @_cmd "sqlite3")
  (#trim! @injection.content 1 1)
  (#set! injection.language "sql")
  (#set! injection.include-children))

(command
  name: (command_name) @_cmd
  argument: (_)
  argument: (string (string_content) @injection.content)
  .
  (#any-basename-of? @_cmd "sqlite3")
  (#set! injection.language "sql")
  (#set! injection.include-children))
