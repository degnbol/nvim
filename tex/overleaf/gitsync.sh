#!/usr/bin/env zsh
# only run for overleaf, make sure we are overleaf
remote=`git remote get-url origin`
[ "${remote[1,25]}" = "https://git.overleaf.com/" ] || return

filename=$1

git add $filename || return 1
git commit -m "local changes: $filename" || return 1
git pull || return 1
git push || return 1

return 0
