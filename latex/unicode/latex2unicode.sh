#!/usr/bin/env zsh
pandoc --filter $0:h/pandoc-unicode-math-from-latex -f latex -t latex
