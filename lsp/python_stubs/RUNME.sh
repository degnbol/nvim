#!/usr/bin/env zsh
cd $0:h
conda create -n python_stubs mypy gemmi
conda activate python_stubs
stubgen --inspect-mode --package gemmi --out .
conda deactivate
conda env remove -n python_stubs
