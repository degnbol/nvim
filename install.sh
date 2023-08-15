#!/usr/bin/env zsh
ln -s ~/dotfiles/config/nvim ~/nvim

$0:h/install-neovim.sh

if [ `uname` = "Darwin" ]; then
    brew install aspell
    ./spell.sh
else
    echo "install aspell then run ./spell.sh"
fi

# https://tree-sitter.github.io/tree-sitter/creating-parsers#installation
cargo install tree-sitter-cli

# ripgrep for telescope to perform searching of words within files
mamba install -yc conda-forge ripgrep pynvim

# consider universal ctags
# brew uninstall ctags
# brew install universal-ctags

# $0:h/install-scala.sh

# for editing jupyter notebooks
if which pipx > /dev/null; then
    pipx install jupytext
else
    mamba install -yc conda-forge jupytext
fi

0:h/tex/unicode/install.sh

# when compiling latex with vimtex, we can use skim instead of default pdf 
# viewer, so that we are at the same location in the document whenever it gets 
# changed. It also auto opens the pdf, and auto updates without having to focus the app first.
brew install skim
echo "Open skim and customize bar to your liking, e.g. hide text and show toggle pane buttons."
# https://dr563105.github.io/blog/skim-vimtex-setup/
echo "Go to sync settings. Enable the reload tickboxes. Put"
echo "Preset -> Custom"
echo "Command -> nvim"
echo "Arguments -> --headless -c \"VimtexInverseSearch %line '%file'\""
echo "Shift+Cmd+click now moves cursor in neovim (if skim was started by neovim)."

# julia LSP
# julia LSP doesn't load info about packages, maybe because it takes too long 
# and a timeout is reached somewhere.
# The hacky solution is to use a precompiled system image as the julia env so everything is already loaded.
# https://discourse.julialang.org/t/neovim-languageserver-jl/37286/83
# first just try to make the regular julia work
$0:h/install.jl

# python LSP "pylsp"
# https://github.com/python-lsp/python-lsp-server
pip install python-lsp-server

# for asciidoc editing.
# tags is to goto definition for xrefs, since there is no LSP, treesitter or similar.
brew install asciidoctor universal-ctags

