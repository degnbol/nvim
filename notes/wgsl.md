# WGSL / Bevy Shader Highlighting

Two tree-sitter parsers installed: `wgsl` (szebniok) and `wgsl_bevy` (tree-sitter-grammars). `wgsl_bevy` is a grammar superset — adds `#import`, `#ifdef`, `#define_import_path`, `virtual`, `as`, `::` on top of base `wgsl`.

**We use the `wgsl` parser** because `wgsl_bevy` cannot parse grouped import syntax (`#import foo::{ bar, baz }`) used by naga_oil in modern Bevy. A single parse error on line 1 cascades and destroys all highlighting. The `wgsl` parser also errors on preprocessor lines but recovers gracefully — standard WGSL highlights correctly.

Bevy preprocessor directives are highlighted via vim regex rules in `syntax/wgsl.vim` as a fallback layer on top of treesitter.

**Upstream issue:** [tree-sitter-grammars/tree-sitter-wgsl-bevy#27](https://github.com/tree-sitter-grammars/tree-sitter-wgsl-bevy/issues/27). If fixed, switch to `wgsl_bevy` via `vim.treesitter.language.register("wgsl_bevy", "wgsl")` and remove `syntax/wgsl.vim`.

**Workaround alternative:** One-import-per-line (`#import bevy_pbr::pbr_fragment::pbr_input_from_standard_material`) makes `wgsl_bevy` parse cleanly.

Note: `wgsl` parser (szebniok) is abandoned (last commit Jan 2023). `wgsl_bevy` is actively maintained.
