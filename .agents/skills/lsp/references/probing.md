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
$ .agents/skills/lsp/hover_probe.py path/to/file.py 106 20 -r path/to/project
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

## Claude Code's LSP tool

The `LSP` tool's hover returns the server's markdown under a position header, so
it shows the server layer. It runs `~/.config/lsp/servers.json`, exported from
this config by `~/.config/lsp/servers.json.sh`. What the export leaves out, and so
where the answers can differ from nvim's, is in the header of
`~/.config/lsp/nvim_servers.lua` and the exclusions in `servers.json.sh` — e.g.
a PEP 723 script gets no per-script client, so its imports resolve from another
environment.
