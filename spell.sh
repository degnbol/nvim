#!/usr/bin/env zsh
# requires aspell. Install on mac with `brew install aspell`
mkdir -p spell
aspell -d en dump master | aspell -l en expand > spell/en.dic
aspell -d da dump master | aspell -l da expand > spell/da.dic

