#!/usr/bin/env zsh
pandoc --filter $0:h/pandoc-unicode-math -f latex -t latex
