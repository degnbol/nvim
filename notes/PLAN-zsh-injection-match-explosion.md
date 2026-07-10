# PLAN: zsh injections.scm — injection silently drops on long commands

## Observed pathology

Miller DSL highlighting (injected into single-quoted `filter`/`put`/`tee`
arguments) works for one `mlr` invocation but not another that is otherwise
equivalent:

```zsh
# works — miller injected
mlr -t --hi sort + uniq -a + filter '$col != ""' ./filename.tsv.gz
# broken — no miller injection
mlr -t --hi --from ./filename.tsv.gz sort + uniq -a + filter '$col != ""'
```

Symptom shape suggested "the `--from` form parses differently" or "the anchor
before the string only matches when a filename trails it". Both are wrong. The
two commands parse to **structurally identical** trees (`command` with a flat
list of `word`/`raw_string` arguments), and the miller injection pattern, run
**in isolation**, matches both. The bug is not in the miller pattern and not in
the parse.

## The real bug

Neovim runs the whole `injections` query for a language through
`Query:iter_matches`, and evaluates injections with tree-sitter's **default
`match_limit = 256`** — the cap on *concurrent in-progress (partial) matches*
during the tree walk. It is hardcoded in core, no override is exposed:

```
runtime/lua/vim/treesitter/languagetree.lua:1125
  self._injection_query:iter_matches(root_node, self._source, start_line, end_line + 1)
```

When the number of simultaneously-open partial matches exceeds 256, tree-sitter
**silently drops matches** — no error, no log. The dropped match happens to
include the miller injection, so highlighting just vanishes.

The count of concurrent partial matches in `queries/zsh/injections.scm`
**grows linearly with the number of arguments in the command being scanned** —
roughly **+24 partials per argument**. So it is not the `--from` form that is
special: *any* command long enough (≈9–11 arguments) that contains an injected
string will silently lose its injection. The two example commands sit on
opposite sides of the 256 threshold (10 vs 11 arguments, plus the string's
position within the command — both feed the partial count).

### Evidence (falsification chain)

| test | result | rules out |
|---|---|---|
| both commands' parse trees | structurally identical | "parses differently" |
| isolated miller pattern (with directives) on broken cmd | 3 captures ✓ | "miller pattern / anchor wrong" |
| full query on broken cmd | 0 captures | — |
| full query, `match_limit` 256 → 1000 | 0 → 3 captures ✓ | confirms it is the limit |
| only **1** completed match at high limit | — | cap is on *concurrent partials*, not completed matches |
| short command, string last | works | "string-last is the trigger" |
| long command, string **not** last | breaks | "string-last is the trigger" (both feed the count) |

### Scaling measurement

Minimum `match_limit` required to keep the injection, vs command length:

| args | current query | `@_interp`→`command_name` only | consolidated + anchored |
|---|---|---|---|
| ~6 | 131 | 47 | 2 |
| ~11 | 275 | 47 | 2 |
| ~15 | 371 | 47 | 2 |
| ~19 | 467 | 47 | 2 |
| ~27 | 659 | 47 | 2 |

The current query is **O(command length)**. A correct injection query is
**O(1)** in the length of unrelated commands — the flat 47 / flat 2 columns.

## Why the current code is flawed

The cost comes entirely from the interpreter patterns (`python -c`, `zsh -c`,
`julia -e`, `Rscript -e`, `node -e`, `lua -e`, plus double-quoted and shell
concatenation variants — ~14 patterns). Two compounding flaws:

1. **Floating (unanchored) capture.** Each interpreter pattern contains
   ```scheme
   [ name: (command_name) @_interp
     argument: (word) @_interp ]      ; matches ANY word in the command
   ```
   The `argument: (word) @_interp` branch (there to support wrappers like
   `uv run python -c`) is unanchored, so the query cursor opens one partial
   match per command argument and keeps it alive until it reaches the string at
   the end. That is O(N) partial matches per pattern, where N = argument count.

2. **Duplication multiplies the constant.** Those ~14 patterns are
   near-identical — they differ *only* in three fields: the interpreter
   basename list, the accepted flag(s), and the target language. Everything else
   (the `[name/argument] @_interp` alternation, the `@_flag . string` anchor,
   the `#trim!`/`#set!` directives) is copy-pasted. So the O(N) floating cost is
   paid ~14 times over → the ~24-per-argument slope.

Two further problems worth calling out:

3. **Silent, non-local failure mode.** Because tree-sitter drops matches with no
   diagnostic, the failure is invisible and looks command-specific. It is
   actually a global property: `jq`, `awk`, `sqlite3` and every other injection
   in this file degrade on long commands too, just with more headroom. Any
   future pattern added to this file lowers the threshold for *all* of them.

4. **`match_limit` is a ceiling, not a fix.** Raising it (even if neovim exposed
   it for injections, which it does not) only moves the cliff. An O(N) query
   will always find a long enough command to break. The limit is what converts
   "query is expensive" into "highlighting silently disappears".

## Fix ideas

The consolidation and the performance fix are the **same change**: remove the
duplication and the floating-capture cost disappears with it. The key move is to
stop capturing the interpreter in the query — resolve it in the Lua directive by
walking the tree, so nothing floats.

### Query: only the flag, anchored to the string (O(1))

Each pattern captures just the command flag adjacent to the injected string:

```scheme
(command
  argument: (word) @_flag
  .
  argument: (raw_string) @injection.content
  (#inject-interp! @_flag)
  (#trim! @injection.content 1 1)
  (#set! injection.include-children))
```

The `word . string` shape is O(1) in concurrent partials — each candidate flag
is tested against its immediate next sibling and discarded at once (the flat
cost the old direct-form flag already paid). There is **no `@_interp` capture**,
so nothing floats, and there is no direct-vs-wrapper split. Patterns reduce to
the string-type variants only: raw_string, double-quoted
`(string (string_content))`, and the `'a'$var'b'` concatenation form — ~4–6
patterns replacing ~14. (The concatenation form is easy to overlook: it is a
real, already-anchored interpreter pattern and must survive the rewrite.)

### Directive: resolve the interpreter by walking back

`#inject-interp!` (registered beside `#inject-by-ext!` in
`plugin/treesitter.lua`, implemented in `lua/utils/treesitter.lua`) does all
interpreter logic in Lua:

1. **Flag-shape gate.** `@_flag` must be a single-dash cluster (`-…`) whose last
   char is a command char — `c` for shells/python, `e` for julia/node/lua/R.
   This subsumes plain `-c`/`-e` *and* combined clusters `-lc`/`-xc`/`-uc`
   (getopt: a clustered command flag must be last, since it consumes the rest of
   the line as the command). A non-flag anchor (`sqlite3 db 'sql'`) or a `--long`
   flag fails here and the match is ignored.
2. **Backward walk.** From `@_flag:prev_sibling()`, skip `-`-prefixed flags
   (`-u`, `-l`, …) and take the first non-flag token — an argument word or the
   `command_name` node.
3. **Table lookup.** `basename → { char, lang }` for python/python3, sh/bash/zsh,
   julia, node, lua, R/Rscript. If the interpreter is in the table and its
   `char` matches the flag's command char, set `injection.language`; otherwise
   leave it unset (same ignore-the-match contract as `#inject-by-ext!`). Adding
   an interpreter is one table line; there are no `--long` flags to list (the
   `--cmd` in the old query was never a real flag on any of these tools).

This runs a handful of node steps per *completed* match in cheap Lua — never a
tree-sitter partial — so it never approaches `match_limit`. It handles every
prefix and wrapper uniformly, because it searches by *what the interpreter is*,
not *what wraps it*:

| command | walk from flag | injects |
|---|---|---|
| `python -u -c 'x'` | `-u` skip → name `python` | python |
| `timeout 180 python3 -u -c 'x'` | `-u` skip → arg `python3` | python |
| `uv run python -c 'x'` | arg `python` | python |
| `env -i bash -c 'x'` | arg `bash` | zsh |
| `zsh -lc 'ls'` | name `zsh`, flag ends `c` | zsh |
| `TZ=UTC uv run python -c` | arg `python` (assignment prefix is a sibling) | python |
| `grep -e 'pat'` | name `grep`, not in table | — |
| `nvim -c 'echo'` | name `nvim`, not in table | — (nvim's own vim pattern fires) |

Concurrent-partial count is **flat at 2 regardless of command length**. This
also recovers the wrapper-with-intervening-flag case (`timeout … python -u -c`)
that a query-only anchored form would have dropped.

### Accepted imperfections

- **Clustered arg-taking middle flag.** `python -mc 'ls'` is `-m c` (module named
  `c`), not `-c`, yet the flag ends in `c` and would misfire. Absurd in practice;
  worst case is cosmetic python highlighting on a non-python string. Not worth
  parsing full getopt semantics in Lua.
- **Interpreter as a bare positional.** `echo hello python -c 'x'` injects python
  (nearest non-flag before `-c` is `python`). This matches the *current* floating
  query's behaviour — no regression. If a real false positive bites, gate the
  walk on the `command_name` being the interpreter itself or in a small wrapper
  allowlist (`timeout`/`time`/`uv`/`env`/…); not worth the maintained list up
  front.

Combined short flags (`zsh -lc`) fall out of the flag-shape gate for free — no
separate follow-up. The old query matched the flag with `#any-of? @_flag "-c"`,
an exact whole-word compare that `-lc` never satisfied; the terminal-char test
in the directive fixes it and combines cleanly with the backward walk.

### Lesser alternatives (documented, not recommended)

- **Keep the interpreter in the query, anchored to the flag.** The version
  grilled before the directive walk: `name:` (direct) + adjacent
  `interp . flag . string` (wrapper) patterns. Truly O(1) too, but drops
  wrapper-with-intervening-flag (`timeout … python -u -c`) and needs a
  direct/wrapper pattern split the walk removes.
- **Collapse `@_interp` to `command_name` only.** Flat 47, but silently drops all
  wrapper support (`uv run python -c` no longer highlights). Simpler diff, worse
  behaviour.
- **Raise `match_limit`.** Not available for injections without patching core,
  and only moves the cliff. Reject.

### Upstream note

Neovim exposes no way to set `match_limit` for injection scanning (hardcoded 256
in `languagetree.lua`). That is the real ceiling; this query is simply the first
in the config to hit it. Worth a local note here and possibly an upstream
report, but the actionable fix is to make the query O(1) so the ceiling is never
approached.

## Scope check

- The miller injection patterns themselves are correct and unchanged.
- The `awk`/`jq`/`sqlite3`/`nvim`/heredoc patterns are not the cause, but they do
  add to the shared partial-match budget; consolidating the interpreter block
  buys headroom for all of them. The generic `word . string` interpreter pattern
  coexists with them: `nvim -c` and `sqlite3 db 'sql'` match it but are rejected
  by the table miss / flag-shape gate, so their own patterns still fire.
- No regression: direct and wrapped interpreter calls, intervening flags, and
  every prefix all keep working, and combined short flags (`zsh -lc`) now work
  too. The only behavioural deltas are the two accepted imperfections above.

## Tests

Extend `tests/plenary/zsh_injections_spec.lua` (helpers `assert_injection` /
`assert_no_injection`; the harness parses via the same 256-limit path
production uses, so long-command explosion reproduces deterministically):

- **Explosion regression** — the original symptom must inject:
  `mlr -t --hi --from ./filename.tsv.gz sort + uniq -a + filter '$col != ""'`
  → miller. Fails on the current query (miller starved), passes after the fix.
  Add a second, longer synthetic command to guard the margin.
- **Combined flags** — `zsh -lc 'ls'` (direct) and `uv run zsh -lc 'ls'`
  (wrapper) inject zsh.
- The existing wrapper / intervening-flag cases (`timeout 180 python3 -u -c`,
  line 160) stay `assert_injection` — the walk keeps them working.
