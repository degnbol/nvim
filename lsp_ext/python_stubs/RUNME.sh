#!/usr/bin/env zsh
# Generate Python type stubs for C-extension packages.
# Stubs are used by basedpyright via the stubPath setting.
cd $0:h

# stub [--skip mod1,mod2] [--static] <package> --with dep [--with ...]
stub() {
  local skip="" mode="--inspect-mode" uv_args=() pkg=""
  while (( $# )); do
    case $1 in
      --skip)   skip=$2; shift 2 ;;
      --static) mode=""; shift ;;
      --with)   uv_args+=(--with $2); shift 2 ;;
      *)        pkg=$1; shift ;;
    esac
  done

  if [[ -n $skip ]]; then
    local mods=($(uv run --no-project $uv_args python3 -c "
import pkgutil, importlib
pkg = importlib.import_module('$pkg')
skip = set('$skip'.split(','))
for _, name, _ in pkgutil.walk_packages(pkg.__path__, '$pkg.'):
    if not any(name == s or name.startswith(s + '.') for s in skip):
        print(name)
"))
    local m_args=(-m $pkg)
    for m in $mods; do m_args+=(-m $m); done
    uv run --no-project --with mypy $uv_args stubgen $mode --out . $m_args
  else
    uv run --no-project --with mypy $uv_args stubgen $mode --out . --package $pkg
  fi
}

stub gemmi    --with gemmi
stub freesasa --with freesasa
stub gudhi    --with gudhi --with scikit-learn --with matplotlib --with pot \
              --skip gudhi.tensorflow
stub Bio      --with biopython --static

# pyrosetta: ~2 GB install, skip if stubs already exist
if [[ ! -d pyrosetta ]]; then
  for pkg in pyrosetta pyrosetta.rosetta pyrosetta.rosetta.std \
             pyrosetta.rosetta.core pyrosetta.rosetta.protocols; do
    stub $pkg --with pyrosetta --with numpy --with dask --with psutil \
              --with billiard --with gitpython --with toolz --with attrs \
              --with distributed
  done
fi
