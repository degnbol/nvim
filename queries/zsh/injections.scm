; extends

; -----------------------------------------------------------------------------
; Quoted-string captures: match the whole (raw_string) / (string) node with
; `#trim!` + include-children, never `(string (string_content))`. A `"…$var…"`
; is parsed by tree-sitter-zsh as string_content + variable_ref siblings, so a
; string_content capture matches once per fragment and injects each as its own
; broken sub-tree; the whole node is one contiguous tree instead. Both quote
; styles delimit with one byte per side, so where double quotes are wanted the
; two node types collapse into `[(raw_string) (string)]` — mlr and awk stay
; raw_string-only on purpose (see their own comments).
;
; `$var` / `${var}` bytes keep their zsh colours: after/queries/zsh/highlights.scm
; re-asserts the base zsh captures over the injected region at priority 101
; (they would otherwise take the injected language's colour). `$(…)` is not
; covered there, so a command substitution inside an injected string takes the
; injected language's colours.
;
; Command-name gates are `#command-is?` (see lua/utils/treesitter.lua) — bar the
; grep patterns, whose predicate resolves the name itself — and it
; resolves wrapper prefixes, so `command jq` / `sudo sqlite3` /
; `env A=1 mlr` count as the wrapped command. tree-sitter-zsh keeps those flat,
; so a wrapper's own tokens become `argument:` siblings of the same `(command)`:
; a pattern that keys off argument *content* is unaffected, but one that counts
; preceding arguments must count from the effective command (`#arg-after?`, used
; by the sqlite3 pattern) rather than from a bare `argument: (_)`.
; -----------------------------------------------------------------------------

; Inject miller DSL into single-quoted strings after put/filter/tee verbs
; in mlr commands. Only raw_string (single-quoted) — double-quoted strings
; have zsh variable expansion which conflicts with miller's $ field refs.
(command
  name: (command_name) @_cmd
  argument: (word) @_verb
  .
  argument: (raw_string) @injection.content
  (#command-is? @_cmd "mlr")
  (#any-of? @_verb "filter" "put" "tee")
  (#trim! @injection.content 1 1)
  (#set! injection.language "miller")
  (#set! injection.include-children))

; Same, but with flags between the verb and the DSL, e.g. `mlr put -q '…'` or
; `mlr put -e '…'`. The dash-flag is anchored adjacent to the string (not just
; "verb somewhere before"): this stops the pattern jumping over a chained verb
; in `mlr put '…' then filter '…'`, where `filter` — not a flag — precedes the
; second string. `-f`/`-s` take a filename / name=value, not DSL, so exclude
; them to avoid injecting a quoted filename as miller.
(command
  name: (command_name) @_cmd
  argument: (word) @_verb
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#command-is? @_cmd "mlr")
  (#any-of? @_verb "filter" "put" "tee")
  (#lua-match? @_flag "^%-")
  (#not-any-of? @_flag "-f" "-s")
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

; Chained verb with flags between it and the DSL, e.g. `… + put -q '…'`.
(command
  name: (command_name) @_verb
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#any-of? @_verb "filter" "put")
  (#lua-match? @_flag "^%-")
  (#not-any-of? @_flag "-f" "-s")
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
; starves the injection out past tree-sitter's hardcoded match_limit (256) on
; long commands.
;
; There is no static injection.language: an off-table interpreter or a
; non-command flag leaves it unset, so the capture is ignored (same contract as
; #inject-by-ext!). That lets `nvim -c` / `sqlite3 db 'sql'` / `grep -e` fall
; through to their own patterns or to no injection.
;
; The `'a'$var'b'` concatenation form needs its own pattern (below).
; -----------------------------------------------------------------------------
(command
  argument: (word) @_flag
  .
  argument: [(raw_string) (string)] @injection.content
  (#inject-interp! @_flag)
  (#trim! @injection.content 1 1)
  (#set! injection.include-children))

; Concatenation form: `zsh -c 'prefix'$var'suffix'`. Each raw_string / string
; fragment is injected independently; variable_ref siblings are highlighted by
; the outer zsh parser. Fragments may not parse as complete code, but most
; highlighting still comes through via error recovery.
(command
  argument: (word) @_flag
  .
  argument: (concatenation
    [(raw_string) (string)] @injection.content)
  (#inject-interp! @_flag)
  (#trim! @injection.content 1 1)
  (#set! injection.include-children))

; Inject awk into the first raw_string argument of awk/gawk/mawk
(command
  name: (command_name) @_cmd
  argument: (raw_string) @injection.content
  (#command-is? @_cmd "awk" "gawk" "mawk")
  (#trim! @injection.content 1 1)
  (#set! injection.language "awk")
  (#set! injection.include-children))

; Inject jq into every raw_string / string argument of jq / gojq. The filter
; is whichever quoted arg(s) the user passes; arg-bearing flags like `--arg`,
; `--argjson` also take quoted values, but those parse as jq string literals
; without breaking, so highlighting them as jq is harmless.
(command
  name: (command_name) @_cmd
  argument: [(raw_string) (string)] @injection.content
  (#command-is? @_cmd "jq" "gojq")
  (#trim! @injection.content 1 1)
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
  argument: [(raw_string) (string)] @injection.content
  (#command-is? @_cmd "nvim" "vim")
  (#any-of? @_flag "-c" "--cmd")
  (#not-lua-match? @injection.content "^['\"]lua\n")
  (#trim! @injection.content 1 1)
  (#set! injection.language "vim")
  (#set! injection.include-children))

; +'...' / +"..." form: concatenation(word("+"), string/raw_string).
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus [(raw_string) (string)] @injection.content)
  (#command-is? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#not-lua-match? @injection.content "^['\"]lua\n")
  (#trim! @injection.content 1 1)
  (#set! injection.language "vim")
  (#set! injection.include-children))

; Multi-line lua: `nvim -c 'lua\nCODE\n'` and `nvim +"lua\nCODE\n"` etc.
; #trim! is a custom directive (see plugin/treesitter.lua) that skips a
; byte-prefix and byte-suffix from the captured node, computing a (row, col,
; byte) range with all three coordinates consistent. `#offset!` can't do this
; because it does naive (row+drow, col+dcol) arithmetic and the col delta
; needed to reach column 0 of the next line varies with the surrounding text.
;
; Skip 5 bytes (`'lua\n` / `"lua\n`), strip the 1-byte closing quote.
; `include-children` keeps the metadata.range from #trim! intact (otherwise
; the range gets masked by the inner string_content/`"` children).
(command
  name: (command_name) @_cmd
  argument: (word) @_flag
  .
  argument: [(raw_string) (string)] @injection.content
  (#command-is? @_cmd "nvim" "vim")
  (#any-of? @_flag "-c" "--cmd")
  (#lua-match? @injection.content "^['\"]lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; +'lua\n...\n' concatenation form
(command
  name: (command_name) @_cmd
  argument: (concatenation (word) @_plus [(raw_string) (string)] @injection.content)
  (#command-is? @_cmd "nvim" "vim")
  (#eq? @_plus "+")
  (#lua-match? @injection.content "^['\"]lua\n")
  (#trim! @injection.content 5 1)
  (#set! injection.language "lua")
  (#set! injection.include-children))

; -----------------------------------------------------------------------------
; Inject regex into the pattern argument of grep / egrep / ggrep / fgrep / rg.
;
; #regex-pattern? (see lua/utils/treesitter.lua) carries the whole decision: it
; resolves the command (so this needs no `#command-is?` duplicating the name
; list, and no @_cmd capture), walks the argument list getopt-style to tell a
; pattern from a glob or a filename, and rejects a pattern whose meaning the
; grammar gets wrong — it models PCRE, so BRE metacharacters and `\d`-style
; escapes read as something else entirely.
;
; A predicate rather than this file's usual `#inject-*!` directive because the
; injected language really is static; a boolean reads better than a directive
; that only ever sets one value.
;
; `grep $'\t…'` is an `ansi_c_string`, a node type the alternation leaves out.
; -----------------------------------------------------------------------------
(command
  argument: [(raw_string) (string)] @injection.content
  (#regex-pattern? @injection.content)
  (#trim! @injection.content 1 1)
  (#set! injection.language "regex")
  (#set! injection.include-children))

; Attached-value form: `grep --regexp='a+'` / `grep -e'a+'` parse as
; concatenation(word, raw_string), so the quoted node is no `argument:` of the
; command. Also catches the fragments of a concatenated pattern (`grep 'a'$x'b'`).
(command
  argument: (concatenation
    [(raw_string) (string)] @injection.content)
  (#regex-pattern? @injection.content)
  (#trim! @injection.content 1 1)
  (#set! injection.language "regex")
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
; `interp <<TAG … TAG` — a command that runs its heredoc body as code because it
; reads its script from stdin (`gnuplot <<GP`, `python <<EOF`, `uv run python -
; <<EOF`). #inject-interp-cmd! takes the whole command node and walks back from
; its last argument (skipping the `-` stdin marker and flags) to the interpreter
; token, then resolves basename → language via the same INTERPRETERS table as
; the `-c`/`-e` flag path (the flag char is irrelevant here); an off-table
; interpreter leaves it unset, so `cat <<EOF` falls through to the
; #inject-by-ext! / heredoc-tag rules.
;
; Two body shapes: a bare command (`python <<EOF`), and a `list` whose last
; command is the interpreter (`cd … && uv run python - <<EOF`).
;
; The heuristic assumes the command's trailing token IS the executed
; interpreter. With no flag to gate on (unlike the `-c`/`-e` path), a command
; whose last argument merely names an interpreter — `echo python <<EOF`,
; `which node <<EOF` — mis-injects its heredoc as that language. Accepted: it's
; inseparable from wrapper support (`timeout 180 python` and `echo python` are
; structurally identical), and such commands feeding a heredoc are rare.
;
; include-children keeps the whole body as one contiguous region: an unquoted
; heredoc (`<<GP`, not `<<'GP'`) expands `$var`, which tree-sitter-zsh splits
; into fragments — capturing them piecewise would inject each as a broken
; sub-tree. The $expansions keep their zsh colours: after/queries/zsh/highlights.scm
; re-asserts the base zsh captures over the injected region at priority 101.
; -----------------------------------------------------------------------------
(redirected_statement
  body: (command) @_cmd
  (heredoc_redirect (heredoc_body) @injection.content)
  (#inject-interp-cmd! @_cmd)
  (#set! injection.include-children))

(redirected_statement
  body: (list (command) @_cmd .)
  (heredoc_redirect (heredoc_body) @injection.content)
  (#inject-interp-cmd! @_cmd)
  (#set! injection.include-children))

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
  (#command-is? @_cmd "nvim" "vim")
  (#eq? @_lflag "-l")
  (#any-of? @_stdin "/dev/stdin" "-")
  (#set! injection.language "lua"))

; -----------------------------------------------------------------------------
; Inject SQL into the query argument of `sqlite3 <db> '<sql>'`.
;
; sqlite3's invocation is `sqlite3 [flags...] <database> <sql>` — the SQL is
; always the LAST argument, preceded by the database path (e.g.
; `"file:zotero.sqlite?immutable=1"`), which is itself frequently quoted. So a
; quoted last argument is the SQL only when the db path precedes it:
; `#arg-after? … 2` requires it to be at least sqlite3's second argument, which
; drops a lone `sqlite3 '<db>'` (interactive) invocation. Counting from the
; *effective* command is what keeps the flag position-independent — a plain
; `argument: (_)` would be filled by a `sudo`/`timeout` wrapper's own tokens, and
; anchoring it adjacent would break on sqlite3's own flags (`-readonly`, `-json`).
; -----------------------------------------------------------------------------
(command
  name: (command_name) @_cmd
  argument: [(raw_string) (string)] @injection.content
  .
  (#command-is? @_cmd "sqlite3")
  (#arg-after? @injection.content @_cmd 2)
  (#trim! @injection.content 1 1)
  (#set! injection.language "sql")
  (#set! injection.include-children))
