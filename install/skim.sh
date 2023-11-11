if [ `uname` = "Darwin" ]; then
    brew install skim
fi
# when compiling latex with vimtex, we can use skim instead of default pdf 
# viewer, so that we are at the same location in the document whenever it gets 
# changed. It also auto opens the pdf, and auto updates without having to focus the app first.
echo "Open skim and customize bar to your liking, e.g. hide text and show toggle pane buttons."
# https://dr563105.github.io/blog/skim-vimtex-setup/
echo "Go to sync settings. Enable the reload tickboxes. Put"
echo "Preset -> Custom"
echo "Command -> nvim"
echo "Arguments -> --headless -c \"VimtexInverseSearch %line '%file'\""
echo "Shift+Cmd+click now moves cursor in neovim (if skim was started by neovim)."

