#!/usr/bin/env zsh
if [ -z "$XDG_CONFIG_HOME" ]; then
    export XDG_CONFIG_HOME=~/dotfiles/config
fi
ln -s $XDG_CONFIG_HOME/nvim ~/nvim
cd $XDG_CONFIG_HOME/nvim/install

./neovim.sh || ./neovim_alt.sh

./spell.sh

# https://tree-sitter.github.io/tree-sitter/creating-parsers#installation
cargo install tree-sitter-cli

# ripgrep for telescope to perform searching of words within files
mamba install -yc conda-forge ripgrep pynvim

# ./scala.sh

# for editing jupyter notebooks
mamba install -yc conda-forge jupytext

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

# for asciidoc editing.
# tags is to goto definition for xrefs, since there is no LSP, treesitter or similar.
brew uninstall ctags
brew install asciidoctor universal-ctags


