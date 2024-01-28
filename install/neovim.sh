#!/usr/bin/env zsh
cd $XDG_CONFIG_HOME/../neovim/
make CMAKE_BUILD_TYPE=RelWithDebInfo
sudo make install
# assuming bin/ and dotfiles/ are in ~
cd ~/bin/
ln -s ../dotfiles/neovim/build/bin/nvim .
