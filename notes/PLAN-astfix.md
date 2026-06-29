# PLAN: `astfix` — ast-grep formatter wrapper + nvim integration

Make the Typst math-symbol ast-grep rulesets usable as a generic "tree-sitter
formatter" both on the CLI and in nvim, keyed by file type, with a notion of a
default rule per language. nvim holds **zero** ast-grep/domain knowledge — it
only shells out to the wrapper on a temp file.

## Verified facts (don't re-litigate)

- `ast-grep scan -c <abs>/sgconfig.yml -r <abs>/rule.yml file.typ` works from
  **any** cwd. `libraryPath: typst.so` in sgconfig resolves relative to the
  **sgconfig.yml's location**, not cwd. → no cwd hack.
- Preview (diff) is the default; `-U` applies in place. → CLI safe by default,
  nvim always passes `-U`.
- `--stdin -U` is **not** a clean stdin→stdout pipe (prints each rule's effect
  against the original separately). → operate on a real file, not stdin.
- ast-grep has no range flag. → range is moot anyway (see nvim § — whole-buffer only).
- `*.so` is a global `.gitignore` rule, so `config/ast-grep/typst.so` is ignored
  at the new path with **no .gitignore change**.
- ast-grep ships an LSP (`ast-grep lsp`) with fix code actions — **not used**: it
  loads a fixed ruleset as persistent diagnostics (can't pick forward/reverse at
  invocation, noisy on every symbol, wants a project-root sgconfig). Shell-out fits.

## Layout (relocation)

Move `config/typst/ast-grep/` → `config/ast-grep/`, rules under a per-extension
subdir. One new top-level folder ever; new languages are subdirs.

```
$XDG_CONFIG_HOME/ast-grep/
├── sgconfig.yml            # registers customLanguages (typst.so), unchanged content
├── typst.so -> …/nvim/site/parser/typst.so   # symlink, gitignored, made at install
├── astfix.zsh             # wrapper  → symlinked bare to ~/.local/bin/astfix
├── _astfix                 # zsh completion → symlinked into zsh/completion/fpath/
├── install.sh              # see below (merges with the existing binary-install one)
├── rules.py                # generator (typst-specific); manual regen only
├── test.sh                 # e2e test, scans rules/typ/* + one wrapper case
└── rules/
    └── typ/                # keyed by file EXTENSION (no ext↔ft map anywhere)
        ├── default.yml     → math-to-unicode.yml   (relative symlink, TRACKED)
        ├── math-to-unicode.yml
        └── math-to-name.yml
```

Decisions locked:
- **Subdir key = extension** (`typ`), not nvim filetype. Both CLI (from filename)
  and nvim (`expand('%:e')`) derive it for free; sgconfig already declares
  `extensions: [typ]`. No mapping table needed.
- **Default = `default.yml` symlink** → `math-to-unicode.yml`. Wrapper needs no
  default logic (resolves the symlink); `--list` skips symlinks so no `default`
  duplicate. The symlink is **relative → portable → tracked** (committed, no
  install step). `rules.py` rewrites the real file in place; symlink survives.
- **`typst.so` symlink is absolute → machine-specific → gitignored**, created at
  install. (Inverse of `default.yml`: that asymmetry is the rule, not an accident.)

## Wrapper: `config/ast-grep/astfix.zsh` (symlinked bare as `astfix`)

zsh script. Resolves its own real dir (`root=${0:A:h}` after symlink resolution)
to find `sgconfig.yml` + `rules/` — works regardless of how it's invoked, no
reliance on `$XDG_CONFIG_HOME`.

Interface:
```
astfix FILE                 # default rule, preview (diff)
astfix -U FILE              # default rule, apply
astfix -r NAME FILE         # named rule (rules/<ext>/NAME.yml)
astfix -r NAME -U FILE      # named rule, apply
astfix --list FILE          # rule basenames for FILE's ext (skip symlinks); no scan
astfix --list               # NO file: union of basenames across all rules/*/ ; for completion
```
Logic:
- `ext = FILE:e`; `dir = $root/rules/$ext`.
- **`--list` branches on argument count**, never erroring:
  - argc 0 → `for f in $root/rules/*/*.yml; [[ -L $f ]] && continue; print ${f:t:r}` | sort -u.
  - else → same over `$dir/*.yml`. **Missing dir → exit 0, empty** (querying an
    unsupported language is benign; lets nvim/zsh completion stay guard-free).
  - nvim passes `expand("%")` (always a file arg, possibly `""`); an empty/unknown
    ext yields empty output, *not* the union — branch on argc, not emptiness.
- **action paths** (`-U`, default scan): missing `dir` or rule file → error +
  non-zero exit. (Acting on an unsupported language is a real error.)
- rule file = `$dir/${NAME:-default}.yml`.
- exec `ast-grep scan -c $root/sgconfig.yml -r <rulefile> [-U] FILE`.
- `set -euo pipefail`; getopts; usage on bad args.
- ponytail: getopts, ~30 lines, no arg-parsing library. Flags precede FILE.

## nvim integration

nvim owns only the temp-file dance — **whole buffer, no range** (visual/partial
selection is out of scope: the rules match `inside: { kind: formula }`, so a
partial slice that splits a `$…$` silently no-ops; whole-buffer is idempotent
and safe).

`lua/astfix.lua` (module):
- `run(rule)`:
  - `ext = vim.fn.expand("%:e")`; if empty → `vim.notify` error, return.
  - tmp = `vim.fn.tempname() .. "." .. ext`; `writefile(nvim_buf_get_lines(0,0,-1), tmp)`.
  - **Build argv without nils** (table literals with embedded `nil` truncate):
    ```lua
    local cmd = { "astfix" }
    if rule then cmd[#cmd+1] = "-r"; cmd[#cmd+1] = rule end
    cmd[#cmd+1] = "-U"; cmd[#cmd+1] = tmp
    ```
  - `vim.system(cmd):wait()`; non-zero `code` → notify `stderr`, return.
  - `nvim_buf_set_lines(0, 0, -1, false, readfile(tmp))`. Single undo step.
- `complete()`: `vim.fn.systemlist({"astfix", "--list", vim.fn.expand("%")})` for
  `:AstFix` arg completion. No guard needed — unknown ext → empty list.

`:AstFix` command (**global**, `plugin/astfix.lua`):
- `nargs = "?"`, `complete = function() return require("astfix").complete() end`.
- callback: `require("astfix").run(args ~= "" and args or nil)`.
- Global is deliberate: `:AstFix` exists everywhere as the discoverable entry; in
  a buffer whose ext has no rules, completion is empty (clear "set up, nothing
  here") and an action errors cleanly. Concrete future consumers: Julia
  inline↔function, LaTeX unicode. `require` is lazy (only in callbacks).

Keymap in `ftplugin/typst.lua` (buffer-local, per-ft):
- `n <leader>tf` → `require("astfix").run()` (whole buffer, default rule).
- No visual keymap. Reverse direction is `:AstFix math-to-name` (rare → command).

## zsh completion: `config/ast-grep/_astfix`

Hand-written `_arguments` (context-dependent rule names → only hand-written gives
this; lives with the tool because it encodes astfix's interface). `#compdef astfix`:
- `-r` → rule names: if a FILE token is already on the line, `astfix --list "$file"`
  (ext-specific); else `astfix --list` (union). Feed to `_describe`/`compadd`.
  (Flags precede FILE, so `-r <TAB>` usually hits the union; narrowing fires on
  re-edits like `-r <TAB> file.typ`. Cheap, harmless.)
- `-U`, `--list` flags; trailing arg `*:file:_files`.
- Verify with `comptest.zsh` per `/completion` skill.

## Install / generator changes

`config/ast-grep/install.sh` — **merge** with the existing one (which installs the
ast-grep binary + regenerates the native `_ast-grep` completion). Ordered:
1. install binary (existing).
2. regenerate native `_ast-grep` (existing).
3. typst.so symlink (inline, ported from the old typst install.sh):
   `[[ -L typst.so ]] && rm typst.so; ln -s "$parser" typst.so`, with the missing-
   parser warning. `# ponytail: typst.so symlink; add a per-language setup.sh only
   if more grammars need parsers.`
4. `ln -s "$PWD/astfix.zsh" ~/.local/bin/astfix`.
5. `ln -s "$PWD/_astfix" "$(git root)/zsh/completion/fpath/_astfix"`.
- **Not run at install**: `rules.py` (output is tracked) and the `default.yml`
  symlink (relative → tracked).

`rules.py`: change `RULES = Path(__file__).parent / "rules"` →
`… / "rules" / "typ"`. (Generator stays typst-specific; manual regen only.)

`test.sh`: `ast-grep scan -r rules/$rs` → `rules/typ/$rs`. **Add one wrapper case**:
`astfix -U` on a `.typ` fixture asserts forward transliteration (the wrapper's
arg-parse / symlink-resolution / ext-dir logic is the new surface needing a check).

## Verification checklist

- [ ] `astfix --list /tmp/x.typ` → `math-to-name`, `math-to-unicode` (no `default`).
- [ ] `astfix --list` (no file) → union across all `rules/*/` (today: same two).
- [ ] `astfix --list /tmp/x.unknownext` → empty, exit 0.
- [ ] `astfix /tmp/x.typ` previews; `astfix -U /tmp/x.typ` applies forward.
- [ ] `astfix -r math-to-name -U /tmp/x.typ` applies reverse.
- [ ] `astfix -U /tmp/x.unknownext` → clean error, non-zero exit.
- [ ] nvim `<leader>tf` transliterates whole buffer; undo is one step.
- [ ] `:AstFix <Tab>` completes `math-to-name`/`math-to-unicode` in a `.typ` buffer;
      empty (no error) in a buffer with no rules for its ext.
- [ ] `_astfix`: `astfix -r <Tab>` offers rule names (`comptest.zsh`).
- [ ] `test.sh` passes after relocation (incl. the wrapper case).

## Out of scope / YAGNI

- No visual/range keymap, no `vim.ui.select` menu — `:AstFix` + completion is the
  menu; whole-buffer is the only safe granularity for formula-enclosed rules. If
  per-formula scoping is ever wanted, the unit is the treesitter `formula` node
  under the cursor, never a line range.
- No buffer-local `:AstFix` per ft — one global command; per-ft is the keymap's job.
- No ft↔ext map, no per-ft nvim config table — extension keying + own-dir
  resolution removes the need.
- No separate `rules.sh` / per-language `setup.sh` — one language today; inline the
  typst.so symlink, split out only when a second grammar needs a parser.
- Generic formatter beyond typst: drop `rules/<ext>/` for a new language; the
  wrapper and nvim helper already handle it. No code change.
