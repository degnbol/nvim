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
nvim +PackerSync

# coc.nvim requires node.js: https://github.com/neoclide/coc.nvim
# if it fails then install to $HOME/bin
curl -sL install-node.vercel.app/lts | bash ||
curl -sL install-node.vercel.app/lts | bash -s -- --prefix=$HOME
# coc needs pynvim
conda install -y pynvim -c conda-forge


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

# when coc is installed then it should automatically install julia support due to the coc-config.json but otherwise run
# :CocInstall coc-julia coc-python coc-sh coc-r-lsp
# NOTE: Set the python intepreter to conda python with command
# :CocCommand python.setInterpreter
# Followed by selecting the miniconda python in the menu that opens.

# for telescope to perform searching of words within files
if $mac; then
    brew install ripgrep
else
    conda install -y -c conda-forge ripgrep
fi

