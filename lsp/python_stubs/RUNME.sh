#!/usr/bin/env zsh
cd $0:h
# isolated env (can also just run the stubgen calls from relevant env
conda create -n python_stubs mypy gemmi biopython freesasa
conda activate python_stubs

# C modules
./stubgen.sh gemmi
./stubgen.sh freesasa
# Bio.PDB.* missing
stubgen --package Bio --out .

# cleanup
conda deactivate
conda env remove -n python_stubs
