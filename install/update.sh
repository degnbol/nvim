#!/usr/bin/env zsh
# update and override miller syntax highlight file
wget -O $XDG_CONFIG_HOME/nvim/syntax/mlr.vim https://raw.githubusercontent.com/johnkerl/miller/main/vim/syntax/mlr.vim
wget -O $XDG_CONFIG_HOME/nvim/lua/surround.lua https://raw.githubusercontent.com/echasnovski/mini.nvim/main/lua/mini/surround.lua

nvim --headless '+Lazy! sync' +qa

# update tree-sitter-cli
command -v rustup > /dev/null && rustup update
cargo install tree-sitter-cli



