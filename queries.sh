#!/usr/bin/env zsh
# after/queries/{language}/highlights.scm seems to have no effect and
# queries/... completely replaces so we extend the default scemas with this

default=~/.local/share/nvim/site/pack/packer/start/nvim-treesitter/queries
ROOT=~/dotfiles/config/nvim/

for after in $ROOT/after/queries/*; do
    lang=$after:t
    mkdir -p "$ROOT/queries/$lang/"
    cat "$default/$lang/highlights.scm" "$after/highlights.scm" > "$ROOT/queries/$lang/highlights.scm"
done

