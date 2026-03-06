---
name: blink
description: blink.cmp completion plugin internals. Use when debugging completion behaviour, modifying keyword handling, trigger logic, or fuzzy matching in blink.cmp.
---

# blink.cmp Internals

## Keyword Character Detection

Blink has **two implementations** for keyword detection, selected at startup:

1. **Rust** (default, via prebuilt binary): Hardcoded regex `[\p{L}0-9_][\p{L}0-9_-]*$` for backward matching, `^[\p{L}0-9_-]+` for forward. Source: `fuzzy/rust/keyword.rs`.
2. **Lua** (fallback): Uses `\k` vim regex but **overrides `iskeyword`** to a fixed `@,48-57,_,-,192-255` before every match, then restores it. Changing `vim.bo.iskeyword` has NO effect. Source: `fuzzy/lua/keyword.lua`.

**Keyword characters in both**: letters, digits, underscore, hyphen. Hyphen cannot start a keyword.

`fuzzy.is_keyword_character(char)` delegates to `get_keyword_range(char, #char, 'prefix')` and checks if start != end. Special case: `-` always returns true.

## Trigger Flow (`completion/trigger/init.lua`)

`on_char_added(char)` checks in order:
1. **Trigger character** (`is_trigger_character`) — resets context (`trigger.context = nil`), creates fresh context, sends new LSP request. Old items discarded.
2. **Keyword character** (`is_keyword_character`) — updates existing context (same context ID if context exists), re-filters existing items.
3. **Neither** — `trigger.hide()`, menu closes.

The trigger character branch fires FIRST. To prevent a character from resetting the context, it must be blocked as a trigger character AND recognised as a keyword character.

### Blocking trigger characters

`show_on_blocked_trigger_characters` (default `{ ' ', '\n', '\t' }`) prevents characters from being treated as trigger characters even if the LSP declares them. Supports functions for filetype-conditional blocking:

```lua
show_on_blocked_trigger_characters = function()
    if vim.bo.filetype == 'typst' then return { ' ', '\n', '\t', ':' } end
    return { ' ', '\n', '\t' }
end,
```

`show_on_x_blocked_trigger_characters` additionally blocks on insert-enter and post-accept (default `{ "'", '"', '(', '{', '[' }`).

### `is_trigger_character` logic

1. Check if char is in `sources.get_trigger_characters()` (from LSP `completionProvider.triggerCharacters`)
2. If char matches `%a` (alphabetic), return false regardless
3. Check blocked lists
4. Return `is_trigger and not is_blocked`

## Context and Bounds

`context.get_bounds(range)` calls `fuzzy.get_keyword_range(line, cursor_col, range)` and converts to 1-indexed `{ line_number, start_col, length }`.

`within_query_bounds(cursor, include_start_bound)` checks if cursor is inside context bounds. Called from `on_cursor_moved`. The `include_start_bound` parameter is true when the character under cursor is a trigger character.

**0-indexed convention**: `get_keyword_range` returns `(start, end)` as 0-indexed byte positions, end exclusive. In Lua, to get the byte at 0-indexed position `pos`, use `line:byte(pos + 1)`. To get the byte just before `start`, use `line:byte(start)` (because Lua 1-indexed `start` = 0-indexed `start - 1`).

## Fuzzy Matching and `guess_keyword_range`

The Rust fuzzy implementation (`fuzzy/rust/fuzzy.rs`) uses `guess_keyword_range` per-item to extend the keyword range backwards through non-keyword characters. This is how path completion works (e.g., `str/tr` matching `str/trim`).

For each completion item, `guess_keyword` computes an extended needle:
- Start with the keyword range from `get_keyword_range`
- Scan backwards from the keyword start
- At each `is_valid_word_boundary` position, check if `line[idx..keyword_start]` matches the start of the item text
- If yes, extend the keyword start backwards

**Valid word boundaries** (`is_valid_word_boundary`): position 0, uppercase after lowercase, non-alphanumeric char, alphabetic after non-alphabetic, digit after non-digit.

This means items with `:` (like `fig:my-figure`) CAN match even when `:` is not a keyword character — the guess logic extends `m` to `fig:m` when it finds `fig:` on the line before the keyword.

## Monkey-patching Keyword Characters

To add keyword characters for specific filetypes (e.g., `:` for typst references), three things are needed:

1. **Block the char as trigger character** — prevents context reset in `on_char_added`
2. **Patch `fuzzy.is_keyword_character`** — makes the char go through the keyword branch
3. **Patch `fuzzy.get_keyword_range`** — extends bounds backwards through the char so `within_query_bounds` stays valid

The patches must run after blink loads. Use a `config` function in the lazy.nvim spec:

```lua
config = function(_, opts)
    require('blink.cmp').setup(opts)
    local fuzzy = require('blink.cmp.fuzzy')
    -- patches here
end,
```

**Current patches** (in `lua/plugins/blink.lua`): `:` is treated as a keyword character in typst buffers, with the keyword range extending backwards through `:` segments (e.g., `fig:my` instead of just `my`).

## LSP Trigger Characters

Check a server's trigger characters at runtime:

```vim
:lua for _, c in ipairs(vim.lsp.get_clients()) do print(c.name, vim.inspect(c.server_capabilities.completionProvider)) end
```

Or headless: `nvim --headless file.typ -c 'sleep 3' -c 'lua ...' -c 'qa' 2>&1`.

Tinymist declares: `#`, `(`, `<`, `,`, `.`, `:`, `/`, `"`, `@`.
