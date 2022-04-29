" reload kitty when we edit the kitty config file
autocmd bufwritepost kitty.conf :silent !~/.config/kitty/reload.sh
