#!/usr/bin/env zsh
# requires aspell. Install on mac with `brew install aspell`
mkdir -p spell
aspell -d en dump master | aspell -l en expand > spell/en.dic
aspell -d da dump master | aspell -l da expand > spell/da.dic
# if `set spelllang=da` or `set spelllang=en,da` then install danish understood by vim spell with
nvim +'set spell|q'
# then follow prompt and install danish utf8

