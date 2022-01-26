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

# I had this error message from which-key: https://github.com/LunarVim/LunarVim/issues/2139
# I fixed it by modifying ~/.local/share/.../start/.../lua/which-key/keys.lua like so: https://github.com/folke/which-key.nvim/pull/231/files
# It is an error that is not patched yet in which-key resulting from breaking changes in nvim 0.7 , so it should eventually be fixed.

# coc.nvim requires node.js: https://github.com/neoclide/coc.nvim
curl -sL install-node.vercel.app/lts | bash
# when coc is installed then you can get julia support with
:CocInstall coc-julia
