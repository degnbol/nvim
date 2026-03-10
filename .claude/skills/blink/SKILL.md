---
name: blink
description: blink.cmp completion plugin internals. Use when debugging completion behaviour, modifying keyword handling, trigger logic, or fuzzy matching in blink.cmp.
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

## References

- [Internals: keyword detection, trigger flow, context/bounds, fuzzy matching, monkey-patching](references/internals.md)
