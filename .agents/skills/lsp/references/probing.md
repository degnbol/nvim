# Seeing what a language server sends

Hover, completion and diagnostics reach the user through two layers: what the
server sends, and what nvim renders. A bug lives in one or the other, and the
float shows only the second — markdown normalised, escapes applied, overloads
concatenated. To split them, ask the server directly.

## hover_probe.py

`../hover_probe.py` speaks LSP to a server over stdio and prints the raw
`contents.value`. It takes the command and settings from
`vim.lsp.config[<server>]`, so `extraPaths`, `stubPath` and whatever
mason-lspconfig merged all apply — no second copy of `lsp/<server>.lua` to drift.

```console
$ .claude/skills/lsp/hover_probe.py path/to/file.py 106 20 -r path/to/project
```

Line and column are 1-based, as an editor displays them. `--json` gives the whole
result, including `contents.kind`. `--server` and `--language` reach any other
`lsp/<name>.lua`.

## Why not drive nvim

`nvim --headless -l` never attaches a client: it skips filetype detection, so
nothing triggers `vim.lsp.enable`. `-c 'luafile …'` with a real file argument does
attach, but then `vim.lsp.buf.hover` gives a float whose buffer holds the rendered
text. `client:request_sync("textDocument/hover", …)` from there returns the raw
payload and is the right call when the question involves nvim's own config
resolution; the probe is faster when it does not.

## cclsp

`mcp__cclsp__get_hover` is an LSP client like nvim: same server and same settings
would mean the same answers. Today it has neither for Python, and only one of the
two is reachable from `~/dotfiles/config/claude/cclsp.json`:

- **The server is configuration.** It runs pylsp (jedi) for `.py`; the `command`
  field can name `basedpyright-langserver --stdio` instead.
- **The settings are not.** cclsp sends `initializationOptions` and nothing else,
  and basedpyright takes `extraPaths`/`stubPath` only from
  `workspace/didChangeConfiguration` or a `pyrightconfig.json`. Passing them as
  `initializationOptions` measures identical to passing nothing.

So a project's `pyrightconfig.json` is what makes every client agree — cclsp,
nvim and the Claude lint hook read it alike. Without one, cclsp resolves imports
from a different environment than nvim, and the symptom is a **blank** hover
rather than an error: a `numpy` attribute in a project with no venv reads as "no
documentation" instead of "wrong environment".
