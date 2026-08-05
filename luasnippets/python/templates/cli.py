"""One line: what the tool produces from what input.

Reads <what data> from stdin (or file args; ±.gz). Note <one-line caveat, if any>.

Writes <format>, one row per <unit>. Columns:
- `col`: what it holds (omit the gloss when the name says it)

Keep this an I/O contract — what is read, what is written. Rationale, provenance, and
design justification are inner thoughts: put them in a comment at the mechanism, or
nowhere. Don't restate what a flag's help text already says.
"""
import argparse
import sys
import csv

from tap import Tap

def f(v: str) -> str:
    """Useful function.

    Docstring is google-style.

    params:
    - v: Interesting input.

    returns: Useful value.
    """
    return v

class Args(Tap):
    cool_value: str = "foo"
    "Name of column with cool value."

    # Short flags are useful
    def configure(self) -> None:
        self.add_argument("-v", "--cool-value")

def main() -> None:
    args = Args(
        underscores_to_dashes=True,
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    ).parse_args()

    # Read TSV
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    if reader.fieldnames is None:
        args.error("No columns found.")
    elif args.cool_value not in reader.fieldnames:
        args.error(f"Column {args.cool_value} not found among {' '.join(reader.fieldnames)}.")

    # Write TSV (DictWriter suits column passthrough; f-string print or csv.writer are fine for fixed scalar columns)
    w = csv.DictWriter(sys.stdout, fieldnames=[*reader.fieldnames, "bar"], delimiter="\t")
    w.writeheader()
    for row in reader:
        bar = f(row[args.cool_value])
        w.writerow({**row, "bar": bar})
