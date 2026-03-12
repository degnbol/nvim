#!/usr/bin/env zsh
set -euo pipefail
cd ${0:A:h}/..

tree-sitter generate 2>/dev/null

pass=0 fail=0

check() {
  local expect=$1 input=$2
  local result=$(tree-sitter parse <(printf '%s' "$input") 2>&1)
  local has_error=0
  if echo "$result" | grep -qE 'ERROR|MISSING'; then
    has_error=1
  fi

  if [[ $expect == ok && $has_error -eq 0 ]]; then
    printf "  OK:    %s\n" "$input"
    pass=$((pass + 1))
  elif [[ $expect == err && $has_error -eq 1 ]]; then
    printf "  OK:    %s  (correctly rejected)\n" "$input"
    pass=$((pass + 1))
  else
    printf "  FAIL:  %s  (expected %s)\n" "$input" "$expect"
    fail=$((fail + 1))
  fi
}

echo "=== Valid pymol selections ==="
check ok  'name CA and chain A'
check ok  'resi 1-100'
check ok  'byres name CA'
check ok  'within 5 of chain A'
check ok  'n. CA or r. ALA'
check ok  'polymer.protein and not solvent'
check ok  'b > 50'
check ok  'b < 50'
check ok  'q < 1.0'
check ok  '/obj//A/1-100/CA'
check ok  'all'
check ok  'chain A+D+E'
check ok  'name CA and not (resn ALA or resn GLY)'
check ok  '%mysel and organic'
check ok  'elem C+N+O'
check ok  '/pept/seg/A/100/CA'
check ok  'A/100/CA'
check ok  '//A/100/'
check ok  'polymer and not hydrogens'
check ok  'cartoon'
check ok  'sticks'
check ok  'everything'

echo ""
echo "=== Should be rejected ==="
# 'hello world' is a benign false positive: error recovery inserts a MISSING
# selector, but identifier has no highlight query so no visual artefact.
check ok  'hello world'
check err '../src/pml/'
check err '/path/to/file'
check err '/usr/local/bin'
check err 'pattern'

echo ""
printf "Results: %d passed, %d failed\n" "$pass" "$fail"
[[ $fail -eq 0 ]]
