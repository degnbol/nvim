# C/C++ language servers: compile flags and missing headers

When goto-definition or hover finds nothing in C/C++, check the compile command
before suspecting the server.

## Where the command comes from

| | clangd | clice |
|---|---|---|
| Database | `compile_commands.json` in the file's parent dirs, or their `build/` | `compile_commands.json` in the workspace root and its immediate subdirs, and dirs above an opened file (when no rule sets `compile_commands`) |
| Static flags | `compile_flags.txt`; `.clangd` `CompileFlags` (static; cannot run pkg-config) | `clice.toml` / `.clice/config.toml` `[[rules]]` (`append`, `compile_commands`, `default_command`) |
| No entry | `clang <file>` plus `initializationOptions.fallbackFlags` | a guessed command (line-1 `inferred-compile-command` diagnostic), or a rule's `default_command` |
| Set from the client | `initializationOptions.compilationDatabasePath` (a dir); per file at runtime: `workspace/didChangeConfiguration` → `settings.compilationDatabaseChanges` | `initializationOptions`; its `rules` *replace* the config file's |

- clice reads its config file and `initializationOptions` once, at start.
  `rules.compile_commands` and `rules.default_command` need a build from
  2026-09-06 or later (clice-io/clice#664).
- `compilationDatabaseChanges` entries must use the realpath: on macOS
  `/tmp/x.c` is `/private/tmp/x.c`, and the other spelling is ignored.
- With only clangd fallback flags, goto-definition into another translation unit
  stops at the header declaration; with a database it reaches the definition.

## Symptoms of a missing include path

- Types from the unresolved header are undefined, so declarations using them are
  invalid and dropped from the AST. Only those symbols fail: goto-definition on
  `static foo_t f(...)`, where `foo_t` comes from the missing header, returns
  nothing, while `static int g(...)` works. `-ferror-limit=0` does not change
  this.
- clang reports only the **first** missing `#include` of a file; later missing
  ones produce no diagnostic. The visible error can name a header unrelated to
  the failing symbol.

## Diagnosing

- `clangd --check=<file>` prints the compile command it chose and every
  diagnostic, without an editor. It does not apply the client's
  `initializationOptions`; pass `--compile-commands-dir=<dir>` for an external
  database.
- clice log: `<cache_dir>/logs/<run>/master.log` (`cache_dir` defaults to
  `<workspace>/.clice`). Its `compile_args:` line shows `source=` per file
  (`Fallback` = guessed).
- Raw replies inside nvim: [probing.md](probing.md).

## This config

Files without a database get flags from `lua/c_fallback_flags.lua` (pkg-config
packages; add a library there), used by `lsp/clangd.lua` (`fallbackFlags`) and
`lsp/clice.lua` (`default_command`, skipped when the project has a clice config).
clice is pinned by the local Mason registry `lua/mason_overrides/`.
