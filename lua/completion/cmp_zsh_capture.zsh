#!/usr/bin/env zsh
# Minimal completion-capture environment for cmp-zsh (blink `zsh` source).
#
# cmp-zsh's own zshrc capture sources the full interactive shell, which here
# means antidote's zle plugins (fzf-tab, syntax-highlighting, autosuggestions)
# and a colored prompt. Those emit escape sequences and perturb init timing,
# racing cmp-zsh's zpty capture sync so terminal control bytes intermittently
# leak as bogus completion items.
#
# Instead: `zsh -f` (no rc, so none of those plugins load) + just the user's
# completion config (fpath, compinit, zstyles). Deterministic, no escapes, and
# custom fpath completions still resolve. $CMP_ZSH_SHARED is cmp-zsh's capture
# hook, resolved from the plugin's runtime path by the nvim config.
zmodload zsh/zpty || { echo 'error: missing module zsh/zpty' >&2; exit 1 }

zpty z zsh -f -i
zpty -w z " source ~/dotfiles/zsh/completion.zsh"
# Drop completion.zsh's matcher-list: each of its match specs is a separate
# completion pass, and cmp-zsh's compadd hook captures every pass, so each
# candidate surfaces once per spec (3× here). blink does its own fuzzy matching,
# so zsh's case-insensitive/substring matcher is redundant anyway.
zpty -w z " zstyle ':completion:*' matcher-list ''"
source ${CMP_ZSH_SHARED:?CMP_ZSH_SHARED not set by the nvim config}
