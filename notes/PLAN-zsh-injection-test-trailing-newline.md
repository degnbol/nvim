# PLAN: zsh_injections_spec — 22 multi-line/heredoc tests fail on trailing newline

## Observed pathology

Running the full spec reports **22 of 90 failures**, all pre-existing and
unrelated to the interpreter match-explosion fix (confirmed: stashing that
change leaves the same 22 red). Every failure is one of two content shapes —
multi-line `-c` strings or heredoc bodies — and every message is the same:

```
no lua injection with text "print(1)". Got: "print(1)\n", "print(1)\n"
```

The three affected groups (22 total):

- `nvim lua` — multi-line `nvim -c 'lua\n…\n'` and the `--cmd` / `+` / double-quoted variants (4)
- `nvim -l heredoc` — `nvim … -l /dev/stdin <<EOF … EOF` and its variants (5)
- `heredoc by file-redirect extension` — `cat > f.<ext> <<EOF … EOF` for every extension, plus the quoted-dest / heredoc-before-redirect / heredoc-tag-as-language cases (13)

Single-line injections (`python3 -c 'print(1)'`, `jq '.a'`, single-line
`nvim -c 'lua print(1)'`) all **pass**.

## The real bug

`assert_injection` compares the injected region's text to the expected string
with **exact equality**:

```lua
for _, entry in ipairs(inj) do
    if entry.text == expected_text then return end
end
```

Multi-line and heredoc injected ranges legitimately **include the trailing
newline** — it is part of the captured node:

- A heredoc body (`(heredoc_body)`) spans up to and including the newline before
  the end tag, so `local x = 1\n` is the correct content.
- A multi-line `-c 'lua\nprint(1)\n'` after `#trim! 5 1` (skip `'lua\n`, strip
  the closing `'`) is `print(1)\n` — the newline is real string content.

So the region text is `"print(1)\n"` but the test expects `"print(1)"`. Single-
line strings have no trailing newline, which is exactly why they pass. **This is
a test-assertion defect, not a query or directive bug.**

### Evidence

Probing the injection trees directly (fresh buffer, exact harness code):

| input | `#trees` | `#regions` | region text |
|---|---|---|---|
| `nvim -c 'lua\nprint(1)\n'` | 1 | 1 | `"print(1)\n"` |
| `cat > /tmp/a.lua <<'EOF'\nlocal x = 1\nEOF` | 1 | 1 | `"local x = 1\n"` |
| `python3 -c 'print(1)'` | 1 | 1 | `"print(1)"` |

Verification of the fix: replacing the compare with
`entry.text:gsub("\n$", "") == expected_text` turns the suite **green (90/0/0)**,
with no other change.

### The doubled `Got: X, X` is a separate, non-blocking artifact

The failure message lists the region twice, but the probe above shows exactly
**one** tree/region per injection. Reproducing `get_injections` in isolation
returns a single entry; the duplicate appears **only in a full-suite run**, so it
is order/state-dependent — parser or buffer state accumulating across the ~90
tests in one nvim instance (`get_injections` deletes its buffer but the tests
share one process, and the config's `FileType` autocmd starts an independent
highlighter parser on each `filetype = "sh.zsh"` assignment).

It does **not** block the primary fix: `assert_injection` already succeeds if
*any* entry matches, and `assert_no_injection` tests currently pass (the
duplicate is same-content, not a cross-test leak). It only adds noise to failure
messages and is a latent hazard for future count-based assertions.

## Fix

### Primary: tolerate a trailing newline in the injection assertion (test-side)

The assertion's intent is "the injected region **is** this code", not "byte-exact
including the line terminator". Strip a single trailing newline before comparing,
in `assert_injection` and in the ad-hoc `vim.tbl_contains(texts, …)` checks
(concatenation fragments, the `nvim -c` multi-string case, the `--arg` jq case)
that do the same exact match inline.

Prefer this over editing ~22 expected strings to embed `\n`: the helper change is
one line, keeps the expectations readable, and does not couple every multi-line
test to the exact newline-inclusion behaviour of the range calculation.

Trade-off: the suite then no longer asserts *where* the injected range ends. If
locking that matters, add **one** targeted test that asserts the trailing newline
is present for a heredoc body, rather than making every test carry it.

### Secondary: investigate the full-suite duplication (harness hygiene)

Lower priority — noise, not incorrectness. Determine whether stale parsers/
highlighters accumulate across tests in the shared process:

- Does disabling the config `FileType` autocmd for the test (or using an unlisted,
  never-current scratch buffer) collapse the duplicate to one entry?
- If so, either tear down the highlighter in `get_injections`, or dedupe the
  collected results by `(lang, text, range)` before returning.

Do not paper over it in `assert_injection`; fix it in `get_injections` so
`assert_no_injection` and any future count assertions stay trustworthy.

## Open question (confirm, non-blocking)

Were these tests ever green, or did a neovim upgrade change injection-range
newline inclusion? The expected strings were written without `\n`, so either
they never passed or the behaviour shifted under them. Worth a quick `git log`
on the spec + nvim changelog check, but it does not change the fix — the trailing
newline is correct content today regardless of when it started being included.

## Scope check

- The query (`queries/zsh/injections.scm`) and directives
  (`#trim!`, `#inject-by-ext!`, `#inject-interp!`) are correct and unchanged —
  the injected ranges are right.
- Only the test harness/spec changes: relax the comparison, optionally dedupe in
  `get_injections`.

## Tests

- After the assertion relaxation, the full spec passes 90/0/0.
- Optionally add one heredoc test that asserts the trailing newline **is**
  present, to pin the range-end behaviour deliberately rather than by accident.
