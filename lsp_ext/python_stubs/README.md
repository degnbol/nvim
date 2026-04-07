Type stubs for C-extension packages that don't ship their own.
Generated with `stubgen --inspect-mode` from mypy via throwaway uv envs.
`basedpyright --createstub` and monkeytype were tried but produced incomplete results.

Run `./RUNME.sh` to regenerate all stubs.

basedpyright finds these via the `stubPath` setting in `lsp/basedpyright.lua`.
Note: `stubPath` is ignored if `[tool.basedpyright]` exists in pyproject.toml
(basedpyright defaults to `typings/` instead). Keep settings in the LSP config.
