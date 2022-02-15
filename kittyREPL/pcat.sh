#!/usr/bin/env zsh
# bracketed paste
cat <(echo -n $'\e[200~') - <(echo $'\e[201~')
