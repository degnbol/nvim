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
    brew install aspell
    ./spell.sh
else
    cd ~/bin
    wget https://github.com/neovim/neovim/releases/download/nightly/nvim-linux64.tar.gz
    tar xzvf nvim-linux64.tar.gz
    ln -s ~/bin/nvim-linux64/bin/nvim ~/bin/nvim
    rm nvim-linux64.tar.gz
    echo "install aspell then run ./spell.sh"
fi

# https://tree-sitter.github.io/tree-sitter/creating-parsers#installation
cargo install tree-sitter-cli


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
brew uninstall ctags
brew install universal-ctags
pipx install virtualenv

# scala
if $mac; then
    brew install coursier/formulas/coursier
    cs setup
    echo 'export PATH="$PATH:$HOME/Library/Application Support/Coursier/bin"' >> ~/.zshrc
else
    curl -fL "https://github.com/VirtusLab/coursier-m1/releases/latest/download/cs-aarch64-pc-linux.gz" | gzip -d > cs
    chmod +x ./cs
    ./cs setup
    echo 'export PATH="$PATH:$HOME/.local/share/coursier/bin"' >> ~/.zshrc
    rm ./cs
fi


# If using coc then after it is installed it should automatically install julia support due to the coc-config.json but otherwise run
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

# for editing jupyter notebooks
pip install jupytext

# when compiling latex with vimtex, we can use skim instead of default pdf 
# viewer, so that we are at the same location in the document whenever it gets 
# changed. It also auto opens the pdf, and auto updates without having to focus the app first.
brew install skim
echo "Open skim and customize bar to your liking, e.g. hide text and show toggle pane buttons."
# https://dr563105.github.io/blog/skim-vimtex-setup/
echo "Go to sync settings and put"
echo "Preset -> Custom"
echo "Command -> nvim"
echo "Arguments -> --headless -c \"VimtexInverseSearch %line '%file'\""
echo "Shift+Cmd+click now moves cursor in neovim."


