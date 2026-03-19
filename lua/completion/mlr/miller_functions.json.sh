#!/usr/bin/env zsh
# Generate miller_functions.json from mlr -F and mlr -K output.
# Usage: ./miller_functions.json.sh > miller_functions.json
set -euo pipefail

python3 -c "
import json, re, subprocess

def run(args):
    r = subprocess.run(args, capture_output=True, text=True)
    return (r.stdout or '') + (r.stderr or '')

# Parse mlr -F: 'name  (class=X #args=N) Description.'
functions = {}
for line in run(['mlr', '-F']).splitlines():
    m = re.match(r'^(\S+)\s+\(class=([a-z_-]+)\s+#args=(\d+)\)\s+(.*)', line)
    if m:
        functions[m.group(1)] = {
            'class': m.group(2),
            'nargs': int(m.group(3)),
            'desc': m.group(4),
        }

# Parse mlr -K: 'keyword: description' (first line of each entry)
keywords = {}
for line in run(['mlr', '-K']).splitlines():
    m = re.match(r'^([a-z_]+):\s+(.*)', line)
    if m:
        keywords[m.group(1)] = m.group(2)

special_variables = {
    'NR': 'Current record number (across all files)',
    'NF': 'Number of fields in current record',
    'FNR': 'Record number within current file',
    'FILENAME': 'Name of current input file',
    'FILENUM': 'Index of current input file (1-based)',
    'M_PI': 'Mathematical constant pi',
    'M_E': 'Mathematical constant e',
    'IPS': 'Input pair separator',
    'IFS': 'Input field separator',
    'IRS': 'Input record separator',
    'OPS': 'Output pair separator',
    'OFS': 'Output field separator',
    'ORS': 'Output record separator',
    'FLATSEP': 'Separator for flattening nested data',
}

data = {
    'functions': dict(sorted(functions.items())),
    'keywords': dict(sorted(keywords.items())),
    'special_variables': dict(sorted(special_variables.items())),
}
print(json.dumps(data, indent=2))
"
