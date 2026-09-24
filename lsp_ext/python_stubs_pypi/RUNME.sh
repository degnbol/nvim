#!/usr/bin/env zsh
# Types for third-party imports in Python files that have no project venv:
# PEP 561 stub packages (scipy-stubs, ...) and typed runtime packages
# (typed-argument-parser). Installed into a uv venv at
# ~/.local/share/python-stubs/, which basedpyright picks up via extraPaths in
# after/lsp/basedpyright.lua.
# Rerun this script to install or upgrade.

set -euo pipefail

VENV=~/.local/share/python-stubs
uv venv --allow-existing "$VENV"
uv pip install --python "$VENV/bin/python" --upgrade \
    scipy-stubs \
    pandas-stubs \
    typed-argument-parser
