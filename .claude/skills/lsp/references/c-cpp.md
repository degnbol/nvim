# C/C++ language servers: compile flags and missing headers

When goto-definition or hover finds nothing in C/C++, check the compile command
before suspecting the server.

## Where the command comes from

| | clangd | clice |
|---|---|---|
| Database | `compile_commands.json` in the file's parent dirs, or their `build/` | `compile_commands.json` in the workspace root and its immediate subdirs, and dirs above an opened file (when no rule sets `compile_commands`) |
| Static flags | `compile_flags.txt`; `.clangd` `CompileFlags` (static; cannot run pkg-config) | `clice.toml` / `.clice/config.toml` `[[rules]]` (`append`, `compile_commands`, `default_command`) |
| No entry | the nearest entry's command; with an empty database, `clang <file>` plus `initializationOptions.fallbackFlags` | a guessed command (line-1 `inferred-compile-command` diagnostic), or a rule's `default_command` |
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

A project without a database of its own (`find_own` in `lua/compile_db.lua`)
gets one generated in `stdpath("cache")/compile_db/<sha256 of root>/`, before
the server starts (the configs' `root_dir`). Include dirs come from the sources'
`#include` lines: project dirs holding a header, then pkg-config flags of the
packages owning headers not found elsewhere. The database is made once per
root and session; `:lsp restart` reuses it. A new `#include` takes effect on
`<localleader>r` (`compile_db.regenerate`), which regenerates and restarts.
Headers no `.pc` file provides need a database or config file in the project.
The own-database check is per opened file: a database two or more levels below
the root is found only for files under it, and the root's generated database can
hide it.

- `after/lsp/clangd.lua`: `compilationDatabasePath` and `fallbackFlags`.
- `after/lsp/clice.lua`: `rules` with `compile_commands` and
  `default_command`, and `cache_dir` in
  `stdpath("cache")/clice/<sha256 of root>/`; all skipped when the project has
  a clice config file.
- clice is pinned by the local Mason registry `lua/mason_overrides/`.
