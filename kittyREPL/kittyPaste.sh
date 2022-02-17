#!/usr/bin/env zsh
# USAGE: pipe STDIN into `./kittyPaste.sh KITTY_WINDOW_ID`
# bracketed paste escape sequences.
# --copy-env means radian will be able to find R home and ipython will be found in conda
cat <(echo -n $'\e[200~') - <(echo $'\e[201~') | kitty @ send-text --stdin --copy-env --match id:$1
