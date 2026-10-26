#!/usr/bin/env zsh
set -euo pipefail
cd ${0:A:h}
nvim -u ../init.lua -l update.lua
./tree-sitter-typc.sh
