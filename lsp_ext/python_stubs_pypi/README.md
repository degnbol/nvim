PEP 561 stub packages from PyPI (`scipy-stubs`, `pandas-stubs`, ...).
Complements `../python_stubs/` (stubgen for C-extension packages without PyPI
stubs — gemmi, gudhi, freesasa, biopython, pyrosetta).

Run `./RUNME.sh` to install or upgrade. Installs into a uv venv at
`~/.local/share/python-stubs/` (not git-tracked). basedpyright finds them via
`extraPaths` in `../../lsp/basedpyright.lua`.

Add a stub package here when:
- the project commonly uses the library (scipy, pandas, ...),
- a community-maintained stub package exists on PyPI (`<name>-stubs`,
  `types-<name>`), and
- it tracks recent versions (check PyPI before adding; outdated stubs are
  worse than no stubs).
