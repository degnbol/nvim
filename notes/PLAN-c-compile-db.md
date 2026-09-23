# Plan: generated compile_commands.json for C/C++ projects without one
Status: reviewed spec

## Background
A C/C++ project with no compilation database, whose configure step fails (e.g. a
`find_package(X REQUIRED)` for an uninstalled library), gives clangd and clice no
include paths. Headers found only through pkg-config (Homebrew
`include/glib-2.0/glib.h`) go unresolved, declarations using their types are
dropped, and goto-definition on those symbols returns nothing. Today
`lua/c_fallback_flags.lua` fixes this with a hand-written package list.

No existing server, generator or plugin infers flags in this case; the
dependency facts that survive a broken configure are the sources' own
`#include` lines. Mapping each unresolved header to the pkg-config package that
owns it replaces the hand-written list.

The output is a real `compile_commands.json` in a cache directory, not fallback
flags: clangd reaches definitions in other translation units only for database
files, both servers take a database directory at initialize
(`compilationDatabasePath`; clice `rules[].compile_commands`), clangd keeps its
index next to that database, and nothing lands in the project. clang reports
only the first missing header per file and clice reads its config only at start,
so flags are computed before the server starts, through the async
`root_dir(bufnr, on_dir)`. Server behaviour: `.claude/skills/lsp/references/c-cpp.md`.

## Changes

### `lua/utils/c_includes.lua`
- `M.includes(lines: string[]) -> { angle: string[], quoted: string[] }` —
  header names of `#include <…>` / `#include "…"` lines, in order, unique
  (Lua patterns over the line-based directive; a directive inside a block
  comment is reported too — harmless unless it names a real header).
- `M.parse_search_dirs(verbose_output: string) -> string[]` — dirs between
  `#include <...> search starts here:` and `End of search list.`, without
  `(framework directory)` entries.
- `M.search_dirs(compiler: string, language: "c"|"c++", callback: fun(dirs: string[]))`
  — `<compiler> -x<language> -E -v -` on empty stdin. Compiler missing → warn
  once, `{}`.
- `M.find_header(header: string, dirs: string[]) -> string|nil`.

### `lua/utils/pkg_config.lua`
- `M.include_dir_index(callback: fun(index: table<string, string[]>))` —
  package → its *own* `-I` dirs: `pkg-config --maximum-traverse-depth=1
  --cflags-only-I <pkg>` (pkgconf flag; Homebrew and Arch ship pkgconf) for
  every `--list-package-names` entry, at most 16 processes at a time. Cached as
  JSON in `stdpath("cache")`, keyed by `PKG_CONFIG_PATH`, `PKG_CONFIG_LIBDIR`
  and the mtimes of the `pc_path` dirs (a `.pc` edited in place is missed).
  pkg-config missing → warn once, `{}`.
- `M.providers(headers: string[], index: table<string, string[]>) -> string[]`
  — for each header, the package whose own dir holds it (ties: longest dir,
  then name); unowned headers skipped. Sorted, unique.
- `M.cflags` unchanged; its `:wait()` is fine inside the async chain (1–3
  packages, a few ms).

### `lua/compile_db.lua`
- `M.find_own(path: string, root: string|nil) -> string|nil` — a database
  either server would find on its own: upward from `path`,
  `compile_commands.json`, `build/compile_commands.json`, `compile_flags.txt`;
  plus `<root>/*/compile_commands.json`. An injected path would hide these.
- `M.cache_dir(root: string) -> string` — `stdpath("cache")/compile_db/` +
  `vim.fn.sha256(root)`.
- `M.collect_files(root: string) -> { sources: string[], headers: string[] }`
  — `vim.fs.dir` walk skipping dot-dirs and `build`, `build-*`, `build_*`;
  sources `c cc cpp cxx c++`, headers `h hh hpp hxx`. Stops after 20000 walked
  entries (warns; returns what it has). Reads in chunks with `vim.schedule`
  between them.
- `M.infer_flags(root: string, files: {sources, headers}, callback: fun(flags: string[]))`
  — 1. project dirs: for each include not found next to its includer, dirs
  under `root` that hold it (`-I`); 2. remaining angle/quoted headers not in
  project dirs ∪ `search_dirs("clang", "c"/"c++")`; 3. `pkg_config.providers`
  over `include_dir_index`; 4. flags = project `-I`s, then
  `pkg_config.cflags(providers)`.
- `M.write(dir: string, root: string, sources: string[], flags: string[])` —
  `dir/compile_commands.json`, one `{ directory = root, file, arguments =
  { "clang", flags…, "-c", file } }` per source (clangd resolves the driver
  in-process; `clang` need not exist).
- `M.generate(root: string, callback: fun(result: { dir: string, flags: string[] }|nil))`
  — composes the three; one run per root at a time (pending callbacks queued);
  result kept in memory per root.
- `M.root_dir(name: string) -> fun(bufnr: integer, on_dir: fun(root: string|nil))`
  — markers from `vim.lsp.config[name].root_markers` at call time;
  `vim.fs.root(bufnr, markers)`; no root → `on_dir(nil)`; `find_own` →
  `on_dir(root)`; a running `name` client for `root` → `on_dir(root)`
  (FileType re-fires do not regenerate; `:LspRestart` does); else `generate`,
  then `on_dir(root)` if the buffer is still valid (also after a failed run).
- `M.result(root: string) -> { dir: string, flags: string[] }|nil` — for
  `before_init`.

### `lsp/clangd.lua`, `lsp/clice.lua`
- `root_dir = function(bufnr, on_dir) require("compile_db").root_dir("<name>")(bufnr, on_dir) end`
  (required inside: these files load at every startup).
- clangd `before_init`: with `result(root)`, set `compilationDatabasePath =
  result.dir` and `fallbackFlags = result.flags` (header-only projects and
  files outside the walk).
- clice `before_init`: no `clice.toml` / `.clice/config.toml` (injected rules
  replace the file's) → `project = { cache_dir = stdpath("cache")/clice }`
  always, and with `result(root)`: `rules = { { compile_commands = { dir },
  default_command = { "clang", flags… } } }`. With a clice config file, clice
  runs on that file alone.

### Delete / update
- Delete `lua/c_fallback_flags.lua` and its wiring.
- `lua/mason_overrides/clice.lua` comment: the pin is for `rules.compile_commands`
  and `rules.default_command`.
- `.claude/skills/lsp/references/c-cpp.md` "This config" section.

### Tests (`tests/plenary/`)
- `c_includes_spec.lua`: `includes` (spacing, angle vs quoted, duplicates);
  `parse_search_dirs` on captured `clang -v` output; `find_header`.
- `pkg_config_spec.lua`: `providers` with a fake index (unowned, tie-break);
  `include_dir_index` gives glib its own dir only, second call is cached
  (pending without glib).
- `compile_db_spec.lua`, on temp projects: `find_own` (upward, `build/`,
  one-level subdir, `compile_flags.txt`); `collect_files` skip rules and cap;
  `infer_flags` puts a project `<config.h>` dir before pkg-config flags and
  maps `<glib.h>` to glib-2.0; `generate` writes one entry per source outside
  the project and runs once for concurrent calls.
- Headless: a two-file project (definition in `a.c`, call in `b.c`, shared
  header including `<glib.h>`) — clangd and clice goto-definition from `b.c`
  reach `a.c`; the project tree gains no files.

## Expected outcome
Opening a C/C++ file in a project without its own database starts clangd and
clice with a generated one: pkg-config libraries it includes resolve without
configuration, cross-file definitions work, and the project tree gets no files.
Projects with their own database (anywhere either server would look) are
unchanged, apart from clice's cache moving out of the tree.

## Non-goals
- Headers no `.pc` file provides.
- Reading build systems (CMake/meson/SCons/autotools) for flags or defines.
- Regenerating while editing; a new `#include` takes effect on `:LspRestart`.
- Include dirs above the server's root; other clangd consumers (cclsp); ObjC
  (`.m` is often MATLAB, and `#import` is not scanned); choosing between
  clangd and clice; pruning old cache entries.

## Open questions
None.
