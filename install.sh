#!/usr/bin/env zsh

# install Paq
git clone --depth=1 https://github.com/savq/paq-nvim.git "${XDG_DATA_HOME:-$HOME/.local/share}"/nvim/site/pack/paqs/start/paq-nvim

# language server currently jedi-language-server for its simplicity and lack of annoying wrong error detection
# installed with pipx
brew install pipx
pipx install jedi-language-server
# then maybe install some of the optional dependencies listed on their git https://github.com/pappasam/jedi-language-server
brew install yarn
yarn global add diagnostic-languageserver

