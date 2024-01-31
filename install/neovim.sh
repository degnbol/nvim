#!/usr/bin/env zsh
cd $XDG_CONFIG_HOME/../neovim/
# install locally to avoid needing sudo and for simpler uninstall according to INSTALL.md
make CMAKE_BUILD_TYPE=RelWithDebInfo CMAKE_EXTRA_FLAGS="-DCMAKE_INSTALL_PREFIX=$HOME/bin/neovim"
make install
# assuming bin/ and dotfiles/ are in ~
cd ~/bin/
ln -s ../dotfiles/neovim/build/bin/nvim .
