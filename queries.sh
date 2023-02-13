#!/usr/bin/env zsh
# after/queries/{language}/highlights.scm seems to have no effect and
# queries/... completely replaces so we extend the default scemas with this

default=~/.local/share/nvim/site/pack/packer/start/nvim-treesitter/queries
ROOT=~/dotfiles/config/nvim/

for after in $ROOT/after/queries/*; do
    lang=$after:t
    mkdir -p "$ROOT/queries/$lang/"
    if [ -f $after/highlights.scm ]; then
        cat {$default/$lang/,$after}/highlights.scm > "$ROOT/queries/$lang/highlights.scm"
    fi
    if [ -f $after/injections.scm ]; then
        cat {$default/$lang/,$after}/injections.scm > "$ROOT/queries/$lang/injections.scm"
    fi
done

