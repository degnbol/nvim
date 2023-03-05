#!/usr/bin/env zsh
# update and override miller syntax highlight file
wget -O syntax/mlr.vim https://raw.githubusercontent.com/johnkerl/miller/main/vim/syntax/mlr.vim
wget -O lua/surround.lua https://raw.githubusercontent.com/echasnovski/mini.nvim/main/lua/mini/surround.lua
rm -f "$XDG_CONFIG_HOME/nvim/plugin/packer_compiled.lua"
# then in nvim :Lazy -> S
