#!/usr/bin/env zsh

# install Paq
# git clone --depth=1 https://github.com/savq/paq-nvim.git "${XDG_DATA_HOME:-$HOME/.local/share}"/nvim/site/pack/paqs/start/paq-nvim
# or install packer
git clone https://github.com/wbthomason/packer.nvim ~/.local/share/nvim/site/pack/packer/start/packer.nvim


# language server currently jedi-language-server for its simplicity and lack of annoying wrong error detection
# installed with pipx
brew install pipx
pipx install jedi-language-server
# then maybe install some of the optional dependencies listed on their git https://github.com/pappasam/jedi-language-server
brew install yarn
yarn global add diagnostic-languageserver

# index files for LSP, see https://github.com/ms-jpq/coq_nvim
brew install universal-ctags
pipx install virtualenv
