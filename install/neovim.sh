#!/usr/bin/env zsh
if [ `uname` = "Darwin" ]; then
    brew install neovim npm # install npm for :LspInstall python that install python support for Lsp. 
    npm install -g neovim # if I do :checkhealth it makes a warning recommending to install this so I did
else
    mkdir -p ~/.local/bin ~/.local/share
    cd ~/.local/share
    wget https://github.com/neovim/neovim/releases/download/nightly/nvim-linux64.tar.gz
    tar xzvf nvim-linux64.tar.gz
    ln -s ~/.local/share/nvim-linux64/bin/nvim ~/.local/bin/nvim
    rm nvim-linux64.tar.gz
fi

