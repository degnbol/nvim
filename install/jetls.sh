#!/usr/bin/env zsh
set -euo pipefail
cd $0:h

# https://github.com/aviatesk/JETLS.jl
# Requires Julia 1.12+ so install that with e.g.
# juliaup add release
# And e.g. switch to it from 1.11 with
# juliaup default release

JETLS=~/.local/share/JETLS.jl
git clone https://github.com/aviatesk/JETLS.jl.git $JETLS
cd $JETLS
julia --project=. -e 'using Pkg; Pkg.instantiate()'

