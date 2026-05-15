#!/usr/bin/env zsh
# PEP 561 stub packages from PyPI (scipy-stubs, pandas-stubs, ...) installed
# into a uv venv at ~/.local/share/python-stubs/. basedpyright picks them up
# via extraPaths in lsp/basedpyright.lua.
# Rerun this script to install or upgrade stubs.

set -euo pipefail

VENV=~/.local/share/python-stubs
uv venv --allow-existing "$VENV"
uv pip install --python "$VENV/bin/python" --upgrade \
    scipy-stubs \
    pandas-stubs
