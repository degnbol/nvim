---
name: blink
description: blink.cmp completion plugin internals and this config's custom completion providers (Miller DSL, PyMOL). Use when debugging completion behaviour, modifying keyword handling, trigger logic, or fuzzy matching in blink.cmp, editing lua/plugins/blink.lua or lua/completion/*, or working on blink.compat/nvim-cmp coexistence.
---

# blink.cmp

Completion plugin with Rust fuzzy matching. Config in `lua/plugins/blink.lua`.

## Key concepts

- **Keyword characters**: letters, digits, underscore, hyphen (hardcoded in Rust). Hyphen cannot start a keyword. Changing `vim.bo.iskeyword` has no effect.
- **Trigger characters**: from LSP servers. Blocked list prevents specific chars from resetting the completion context.
- **Fuzzy matching**: `guess_keyword_range` extends backwards through non-keyword chars (enables path completion like `str/tr` matching `str/trim`).

## LSP trigger characters

Check at runtime:
```vim
:lua for _, c in ipairs(vim.lsp.get_clients()) do print(c.name, vim.inspect(c.server_capabilities.completionProvider)) end
```

Tinymist declares: `#`, `(`, `<`, `,`, `.`, `:`, `/`, `"`, `@`.

## Snippet integration

LuaSnip snippets appear via the built-in `luasnip` source. Snippet jump keymaps (`<C-.>`/`<C-,>`) are shared between blink and LuaSnip. See the [snippets skill](.claude/skills/snippets/SKILL.md) for node types, custom utilities, and the add-snippet workflow.

**Buffer-local keymaps caveat:** Blink applies its keymaps (including snippet jump in `i`/`s` modes) as buffer-local mappings on `InsertEnter`. When expanding snippets programmatically via `ls.snip_expand` on a newly opened buffer (without entering insert mode first), these keymaps won't exist. Fix: fire `vim.api.nvim_exec_autocmds('InsertEnter', { buffer = 0 })` before expanding.

## Community source tips

- **Comma-separated values**: `guess_keyword_range` extends backwards through commas, filtering out items that don't match the extended text. Fix by returning `","` from `get_trigger_characters()` — blink resets context after commas so the keyword starts fresh. The source must self-gate (return empty when not in context) since the trigger fires globally. See [internals](references/internals.md) for details.

- **Dash-prefixed items** (CLI flags like `-f`, `--output`): hyphen can't start a keyword, so typing `-f` gives keyword `f`. Accepting item `-f` replaces only `f` → `--f`. Fix: use `textEdit` with a range that includes the typed dash prefix. Compute by matching `%-[-a-zA-Z]*$` on text before cursor.

## blink.compat / nvim-cmp coexistence

blink.compat with `impersonate_nvim_cmp = true` replaces
`package.loaded["cmp"]` with a mock. `require("cmp")` returns the mock even
after clearing package.loaded (blink.compat has its own `lua/cmp/init.lua`). In
`cmp.lua`, check `package.loaded["blink.cmp"]` at runtime and skip cmp setup
when blink is active.

## Config completion providers

Custom blink.cmp providers built for this config. Data files are generated, not
hand-maintained — regenerate after the upstream tool upgrades.

### Miller DSL (`blink_mlr`)

Inside `put`/`filter` DSL strings:
- After `$` → column names from referenced input files.
- Otherwise → DSL builtin functions (223), keywords (40), special variables (14).

Data generated from `mlr -F` and `mlr -K` by
`lua/completion/mlr/miller_functions.json.sh` → `miller_functions.json`.
Regenerate after Miller upgrades.

### PyMOL

Three independent data sources (not derivations of each other):

1. **`lsp_ext/pymol-open-source/modules/pymol/`** — vendored pymol source, added
   to `$PYTHONPATH` in `plugin/paths.lua` so basedpyright resolves `from pymol
   import cmd` (type info, hover, go-to-definition).
2. **`lua/completion/pymol/pymol.html`** — command reference from
   `cmd.write_html_ref()` (usage/args/examples). Not currently wired into
   completions.
3. **`lua/completion/pymol/pymol_settings_descriptions/*.md`** — scraped from the
   PyMOL Wiki. Used by the `pymol_settings` provider
   (`lua/completion/pymol/blink_pymol_settings.lua`), which completes settings
   inside `set()`/`cmd.set()` calls.

`lua/completion/pymol/blink_pymol_select.lua` completes selection keywords,
builtins, operators, and representation names inside Python strings. Activates
only where the `pymol_select` treesitter injection is active (checks
`parser:language_for_range():lang()`) — no completions in docstrings or
non-pymol strings.

Loading is conditional: `ftplugin/python.lua` scans the first 10 lines for
`import.*pymol` (or `<localleader>+` manually), setting `vim.g.loaded_pymol`
which gates the blink providers.

## References

- [Internals: keyword detection, trigger flow, context/bounds, fuzzy matching, monkey-patching](references/internals.md)
