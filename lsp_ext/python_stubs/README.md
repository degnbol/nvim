Type stubs for C-extension packages that don't ship their own.
Generated with `stubgen --inspect-mode` from mypy via throwaway uv envs.
`basedpyright --createstub` and monkeytype were tried but produced incomplete results.
`RUNME.sh` then runs every tree but pyrosetta through `../stub_patches/patch_stubs.py`, which
writes the runtime docstrings in (stubgen writes none). See
`../stub_patches/README.md`.

RDKit is a Boost.Python binary that ships its own typed stubs bundled as a
`rdkit-stubs/` sibling of the installed package (`make stubs`). Those have real
function signatures but five defects, since pybind11-stubgen can't parse
Boost.Python's own signatures and doesn't record what rdkit's Python layer does
to itself at import: raw C++ in docstrings, enum members named after Python
keywords emitted as invalid syntax, untyped `(*args, **kwargs)` properties,
module-level names that only exist at runtime, and methods patched in Python
over a C++ original. `RUNME.sh` copies the stubs here, and the patch sequence
repairs those with `../stub_patches/patch_pybind_stubs.py`.

The modules rdkit's stub build left out entirely — the ones whose bodies sit
inside a `try:` import of a dependency that build lacked (pandas for
`PandasTools`, IPython for `Draw.IPythonConsole`) — get a stub from `stubgen` in
`RUNME.sh` first, which the bundle is then merged over. pyright does not
fall back to source for a module missing from a stub package, so without them
`from rdkit.Chem import PandasTools` is an unknown import symbol.

Vendored as a *complete* package so this `stubPath` copy wins over any
`rdkit-stubs` installed in a project's environment (an installed one otherwise
takes precedence and shows the unrepaired C++/untyped forms).

Run `./RUNME.sh` to regenerate all stubs.

basedpyright finds these via the `stubPath` setting in `after/lsp/basedpyright.lua`.
`stubPath` is ignored in a project owning a `pyrightconfig.json` or a
`[tool.basedpyright]` section. Keep settings in the LSP config.

For libraries with community-maintained PEP 561 stubs on PyPI (scipy-stubs,
pandas-stubs, ...), use `../python_stubs_pypi/` instead — those are curated
and richer than stubgen output.
