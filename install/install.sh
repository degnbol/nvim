#!/usr/bin/env zsh
set -euo pipefail
cd $0:h/..
# If detached HEAD
git switch main
git submodule update --init .
cd install/

# Installing from source to get :restart and vim.pack
./neovim_source.sh

# We are tracking the dic files in git so no need to build them.
# ./spell.sh

command -v npm > /dev/null || echo "Install npm to get LSP installs"

# tree-sitter CLI required for parser compilation (generate + build)
# https://tree-sitter.github.io/tree-sitter/creating-parsers#installation
if ! command -v tree-sitter > /dev/null; then
    command -v cargo > /dev/null || ~/dotfiles/config/cargo/install.sh
    cargo install --locked tree-sitter-cli
fi
# Install treesitter parsers (waits for async compilation to finish)
nvim --headless +"lua require('nvim-treesitter')._install_task:wait()" +qa

# ripgrep for telescope to perform searching of words within files
# conda install -yc conda-forge ripgrep pynvim
# uv tool install ripgrep pynvim

# ./scala.sh

# for editing jupyter notebooks
# conda install -yc conda-forge jupytext

../tex/unicode/install.sh

# tectonic: self-contained XeTeX engine snacks.image prefers first for inline
# math rendering (lua/plugins/pickers.lua). XeTeX reads unicode math natively.
if ! command -v tectonic > /dev/null; then
    source ~/dotfiles/shell/inst.sh   # inst() → brew / pacman / …
    inst tectonic
fi

# Pre-warm tectonic's package bundle. The first compile that loads unicode-math
# (the snacks inline-math template) downloads fontspec/unicode-math/
# latinmodern-math.otf — tens of seconds — so the first math buffer opened on a
# fresh machine looks like it fails to render. Fetch it once, here. Mirrors the
# template's package set (xcolor + lua/plugins/pickers.lua, sorted) and snacks'
# convert invocation (convert.lua). Sentinel makes it a one-time op.
prewarm=~/.cache/nvim/tectonic-prewarmed
if command -v tectonic > /dev/null && [[ ! -f $prewarm ]]; then
    tmp=$(mktemp -d)
    cat > "$tmp/warm.tex" <<'EOF'
\documentclass[preview,border=0pt,varwidth,12pt]{standalone}
\usepackage{amscd, amsmath, mathtools, unicode-math, xcolor}
\begin{document}$α$\end{document}
EOF
    if tectonic -Z continue-on-errors --outdir "$tmp" "$tmp/warm.tex" > /dev/null 2>&1; then
        mkdir -p $prewarm:h && touch $prewarm
    fi
    rm -r "$tmp"
fi

# Blender Python stubs for LSP completion in .blend.py files
uv pip install --target ../lsp_ext/blender-stubs fake-bpy-module-latest

# julia LSP
# julia LSP doesn't load info about packages, maybe because it takes too long 
# and a timeout is reached somewhere.
# The hacky solution is to use a precompiled system image as the julia env so everything is already loaded.
# https://discourse.julialang.org/t/neovim-languageserver-jl/37286/83
# first just try to make the regular julia work
./install.jl

# python LSP "pylsp"
# https://github.com/python-lsp/python-lsp-server
pip install python-lsp-server

# As fallback to LSP for goto def using CTRL-].
# E.g. for asciidoc.
# tags is to goto definition for xrefs, since there is no LSP, treesitter or similar.
if command -v brew > /dev/null; then
    brew uninstall ctags
    brew install asciidoctor universal-ctags
elif command -v pacman > /dev/null; then
    sudo pacman -Sy --no-confirm universal-ctags asciidoctor asciidoctor-pdf
else
    echo "Install universal-ctags"
    return 1
fi
