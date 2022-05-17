See `:h initialization` for order of calls under startup.
Broadly
 - init.lua runs lua/ files.
 - ftplugin/
 - plugins has a similar structure so may run code at all kinds of timings.
   See their code in `~/.local/share/nvim/site/pack/packer/start/`
 - after/
 - after/syntax/
