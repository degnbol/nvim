#!/usr/bin/env zsh
# update and override miller syntax highlight file
wget -O syntax/mlr.vim https://raw.githubusercontent.com/johnkerl/miller/main/vim/syntax/mlr.vim

wget -O lua/surround.lua https://raw.githubusercontent.com/echasnovski/mini.nvim/main/lua/mini/surround.lua

$0:h/PackerSync.sh

