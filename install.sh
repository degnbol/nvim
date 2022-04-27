#!/usr/bin/env zsh
if [ `uname` = "Darwin" ]; then
    mac=true
else
    mac=false
fi

if $mac; then
    # --HEAD for development version instead of stable necessary for some lua config
    brew install --HEAD luajit
    brew install --HEAD neovim
    brew install npm # install npm for :LspInstall python that install python support for Lsp. 
    npm install -g neovim # if I do :checkhealth it makes a warning recommending to install this so I did
else
    cd ~/bin
    wget https://github.com/neovim/neovim/releases/download/nightly/nvim-linux64.tar.gz
    tar xzvf nvim-linux64.tar.gz
    ln -s ~/bin/nvim-linux64/bin/nvim ~/bin/nvim
    rm nvim-linux64.tar.gz
fi

# Install packer
git clone https://github.com/wbthomason/packer.nvim ~/.local/share/nvim/site/pack/packer/start/packer.nvim

# coc.nvim requires node.js: https://github.com/neoclide/coc.nvim
curl -sL install-node.vercel.app/lts | bash || # if it fails then install to $HOME/bin
curl -sL install-node.vercel.app/lts | bash -s -- --prefix=$HOME

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

# when coc is installed then it should automatically install julia support due to the coc-config.json but otherwise run
:CocInstall coc-julia coc-python coc-sh coc-r-lsp

# for telescope to perform searching of words within files
if $mac; then
    brew install ripgrep
else
    conda install -c conda-forge ripgrep
fi
