#!/usr/bin/env zsh
mkdir -p ~/.local/share/
cd ~/.local/share/

git clone --depth=1 https://github.com/neovim/neovim || return 1
cd neovim

make CMAKE_BUILD_TYPE=RelWithDebInfo
sudo make install
