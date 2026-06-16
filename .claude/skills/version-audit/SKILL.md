---
name: version-audit
description: Neovim version audit procedure for this config. Use when neovim has been updated, or when `nvim --version` is newer than the "Last audited" version below — scans the config for deprecated/removed APIs and new built-in replacements, and keeps the nvim-api-guard hook's denylist in sync.
---

# Neovim Version Audit

Last audited: 0.12.1

## Procedure

Run this audit when the user says neovim has been updated, or when
`nvim --version` shows a version newer than "Last audited" above.

### 1. Identify what changed

```
$VIMRUNTIME = output of: nvim --headless -c 'lua print(vim.env.VIMRUNTIME)' -c 'qa' 2>&1
```

Read these files for every version between "Last audited" and current:

| Source | Contains |
|--------|----------|
| `$VIMRUNTIME/doc/news.txt` | Current dev cycle (breaking, features, removed) |
| `$VIMRUNTIME/doc/news-0.XX.txt` | Past release cycles (one file per major) |
| `$VIMRUNTIME/doc/deprecated.txt` | All deprecations, tagged `deprecated-0.XX` |

**Cross-version jumps**: if last audited is 0.11.x and current is 0.13.x, read
`news-0.12.txt`, `news-0.13.txt` (if exists), `news.txt`, and the
`deprecated-0.12` + `deprecated-0.13` sections of `deprecated.txt`.

**Same major, different patch** (0.12.0 → 0.12.1): `news.txt` covers the
entire 0.12 cycle with no per-patch breakdown. Re-reading is idempotent.

### 2. Scan the config

Grep `~/nvim/` (excluding `modules/lazy_repro/`) for:

- APIs listed in `BREAKING CHANGES` and `REMOVED FEATURES`
- APIs listed in the relevant `deprecated-0.XX` sections
- `vim.fn.*` calls that now have Lua-native replacements (from `NEW FEATURES`)
- Removed or renamed options, events, highlight groups, commands
- Deprecated treesitter APIs, LSP methods, or diagnostic interfaces
- New built-in features that replace patterns currently done manually in the config

Report each finding with: file, line, what to change, and which news section
mentions it.

### 3. Update this file

- Update "Last audited" to the current version.
- Update the tables below (add new entries, remove any that were wrong).
- Note any breaking changes that required config fixes in the audit log.

### 4. Update hooks

Glob `~/.config/claude/hooks/nvim*.sh` for all neovim-related hooks. Read each
one to understand what data it relies on from this file, then update both the
tables here and the hook's hardcoded lists to stay in sync.

## vim.fn calls with Lua-native replacements

Used by the `nvim-api-guard` hook to warn on vim.fn.* calls that have
replacements. Only includes functions where the replacement is stable and
complete (not partial or experimental).

| vim.fn | Replacement | Since |
|--------|-------------|-------|
| `bufnr()` | `vim.api.nvim_get_current_buf()` | 0.5 |
| `winbufnr(win)` | `vim.api.nvim_win_get_buf(win)` | 0.5 |
| `bufname(buf)` | `vim.api.nvim_buf_get_name(buf)` | 0.5 |
| `winnr()` | `vim.api.nvim_get_current_win()` | 0.5 |
| `tabpagenr()` | `vim.api.nvim_get_current_tabpage()` | 0.5 |
| `line(".")` | `vim.api.nvim_win_get_cursor(0)[1]` | 0.5 |
| `col(".")` | `vim.api.nvim_win_get_cursor(0)[2] + 1` | 0.5 |
| `getline(n)` | `vim.api.nvim_buf_get_lines(buf, n-1, n, true)[1]` | 0.5 |
| `setline(n, s)` | `vim.api.nvim_buf_set_lines(buf, n-1, n, true, {s})` | 0.5 |
| `expand("%:p")` | `vim.api.nvim_buf_get_name(0)` | 0.5 |
| `json_encode(v)` | `vim.json.encode(v)` | 0.6 |
| `json_decode(s)` | `vim.json.decode(s)` | 0.6 |
| `delete(path)` | `vim.uv.fs_unlink(path)` | 0.10 |
| `mkdir(path, "p")` | `vim.fn.mkdir()` still simpler for "p" flag | — |
| `filereadable(p)` | `vim.uv.fs_stat(p)` | 0.10 |
| `isdirectory(p)` | `vim.uv.fs_stat(p)` check `.type == "directory"` | 0.10 |
| `system(cmd)` | `vim.system(cmd):wait()` | 0.10 |
| `systemlist(cmd)` | `vim.split(vim.system(cmd):wait().stdout, '\n', { trimempty = true })` | 0.10 |
| `jobstart(cmd)` | `vim.system(cmd)` or `vim.uv.spawn()` | 0.10 |
| `timer_start(ms, fn)` | `vim.defer_fn(fn, ms)` or `vim.uv.new_timer()` | 0.10 |
| `glob(pat)` | `vim.fs.glob(pat)` | 0.11 |

### Functions with NO replacement (do not warn)

These are vimscript-only functions with no `vim.api` or `vim.*` equivalent.
The hook should ignore them.

`stdpath`, `has`, `confirm`, `input`, `feedkeys`, `getchar`, `getcharstr`,
`getreg`, `setreg`, `getregtype`, `expand` (general), `fnamemodify`,
`matchadd`, `matchdelete`, `search`, `searchpos`, `searchpair`,
`searchpairpos`, `synID`, `synIDattr`, `synIDtrans`, `pumvisible`,
`complete`, `complete_info`, `mode`, `visualmode`, `getpos`, `setpos`,
`getcharpos`, `setcharpos`, `cursor`, `winsaveview`, `winrestview`,
`getwininfo`, `getbufinfo`, `gettabinfo`, `argc`, `argv`, `tempname`,
`shellescape`, `fnameescape`, `escape`, `nr2char`, `char2nr`, `strdisplaywidth`,
`strchars`, `byteidx`, `charidx`, `readfile`, `writefile`, `rename`,
`getftype`, `resolve`, `simplify`, `pathshorten`, `executable`,
`exepath`, `environ`, `getenv`, `setenv`.

## Audit log

### 0.12.1 (2026-04-13) — initial audit

No config changes needed. New deprecation: `nvim_set_decoration_provider`
`on_line` → use `on_range` instead (not used in config).
