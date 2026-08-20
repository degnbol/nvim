Type stubs for C-extension packages that don't ship their own.
Generated with `stubgen --inspect-mode` from mypy via throwaway uv envs.
`basedpyright --createstub` and monkeytype were tried but produced incomplete results.

RDKit is a Boost.Python binary that ships its own typed stubs bundled as a
`rdkit-stubs/` sibling of the installed package (`make stubs`). Those have real
function signatures but five defects, since pybind11-stubgen can't parse
Boost.Python's own signatures and doesn't record what rdkit's Python layer does
to itself at import: raw C++ in docstrings, enum members named after Python
keywords emitted as invalid syntax, untyped `(*args, **kwargs)` properties,
module-level names that only exist at runtime, and methods patched in Python
over a C++ original. `fix_pybind_stubs.py` copies the stubs here and fixes those
— see its docstring, with `test_fix_pybind_stubs.py` covering the transforms.
Each rewritten file opens with a `# fix_pybind_stubs: <package> <version> <sha>`
comment, recording what produced it without reaching pyright's hover. `<sha>` is
the first 8 hex digits of the script's own sha256, so a consumer can tell a tree
patched by an older version of it from a current one.

The modules rdkit's stub build left out entirely — the ones whose bodies sit
inside a `try:` import of a dependency that build lacked (pandas for
`PandasTools`, IPython for `Draw.IPythonConsole`) — get a stub from `stubgen` in
`RUNME.sh` first, which the fixer then merges the bundle over. pyright does not
fall back to source for a module missing from a stub package, so without them
`from rdkit.Chem import PandasTools` is an unknown import symbol.

Vendored as a *complete* package so this `stubPath` copy wins over any
`rdkit-stubs` installed in a project's environment (an installed one otherwise
takes precedence and shows the unfixed C++/untyped forms).

Run `./RUNME.sh` to regenerate all stubs.

basedpyright finds these via the `stubPath` setting in `lsp/basedpyright.lua`.
Note: `stubPath` is ignored if `[tool.basedpyright]` exists in pyproject.toml
(basedpyright defaults to `typings/` instead). Keep settings in the LSP config.
Such projects only ever read the environment's own `rdkit-stubs/`, which
`fix_pybind_stubs.py --in-place` fixes where it is installed (atomic
whole-tree swap, so a reader never sees a half-rewritten tree).

For libraries with community-maintained PEP 561 stubs on PyPI (scipy-stubs,
pandas-stubs, ...), use `../python_stubs_pypi/` instead — those are curated
and richer than stubgen output.
