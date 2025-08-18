#!/usr/bin/env zsh
cd $0:h

# Installing from source to get :restart and vim.pack
./neovim_source.sh

# We are tracking the dic files in git so no need to build them.
# ./spell.sh

# https://tree-sitter.github.io/tree-sitter/creating-parsers#installation
cargo install tree-sitter-cli

# ripgrep for telescope to perform searching of words within files
conda install -yc conda-forge ripgrep pynvim

# ./scala.sh

# for editing jupyter notebooks
conda install -yc conda-forge jupytext

../tex/unicode/install.sh

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
if [ `uname` = "Darwin" ]; then
    brew uninstall ctags
    brew install asciidoctor universal-ctags
elif command -v pacman > /dev/null; then
    sudo pacman -Sy universal-ctags
else
    echo "Install universal-ctags"
    return 1
fi
