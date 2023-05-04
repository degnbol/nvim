See `:h initialization` for order of calls under startup.
Broadly
 - init.lua runs lua/ files.
 - ftplugin/
 - plugins has a similar structure so may run code at all kinds of timings.
   See their code in `~/.local/share/nvim/site/pack/packer/start/`
 - seems after/ is never run?
 - after/.../

ctrl+f(oward) and ctrl+b(ack) to navigate hover windows.
K to jump to definition in help docs and ctrl+o to go back. With lsp attached 
we remap K to show hover, and double K to move cursor into hover.
ctrl+e(extra) and ctrl+y(ester) to move screen.
zz to center line in screen, zt to put it at top of screen and zb for bottom.

miller (mlr) syntax highlighting added from https://github.com/johnkerl/miller/tree/main/vim
to ftdetect/mlr.vim and syntax/mlr.vim
It can also be done as a plugin with packer:
use {"johnkerl/miller", ft="mlr", rtp="vim"}
but it seemed the ftdetect was overridden to "conf" and it requires all of 
miller git be cloned so a bit much. ftdetect/mlr.vim probs won't need updating, 
other files probs won't be added and I don't edit miller files that much so it 
seems fine to update manually with copy paste with
./update.sh

word count: g then ctrl+g

multiline search with :S in-place of / (see plugin/search.vim)

Last macro can be repeated with @@
Last command with @: which also can take a count.

vim has builtin completion of a whole line.
Start to write a line written somewhere before in the file.
Press <C-x><C-l> (:h i_CTRL-X_CTRL-L) and you get completion window.

ctrl+] is default vim go-to-definition.

