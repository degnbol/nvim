#!/usr/bin/env zsh
# Build the SeniorMars typst code-mode parser (typc) from source into neovim's
# site dir, for ```typc ...``` injection.
set -euo pipefail

SRC="${XDG_DATA_HOME:-$HOME/.local/share}/tree-sitter-typc-src"
SITE="${XDG_DATA_HOME:-$HOME/.local/share}/nvim/site"

command -v tree-sitter > /dev/null || { echo "tree-sitter CLI required (install/install.sh)" >&2; exit 1; }
command -v node > /dev/null        || { echo "node required for tree-sitter generate" >&2; exit 1; }

if [ -d "$SRC/.git" ]; then
    git -C "$SRC" fetch -q --depth 1 origin HEAD
    git -C "$SRC" reset -q --hard FETCH_HEAD
else
    git clone -q --depth 1 https://github.com/SeniorMars/tree-sitter-typst "$SRC"
fi
cd "$SRC"

# Generate only the code variant (markup/math variants unused).
# Builds to build/typc.
sh scripts/generate-variants.sh code

# Install the artifact: core loads the parser by name from site/parser/.
mkdir -p "$SITE/parser" "$SITE/queries/typc"
cc -shared -fPIC -O2 -I build/typc/src -o "$SITE/parser/typc.so" \
    build/typc/src/parser.c build/typc/src/scanner.c
cp build/typc/queries/typst/*.scm "$SITE/queries/typc/"
