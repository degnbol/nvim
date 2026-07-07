"""Title.

Usage details that aren't better described by a flag's help text.
"""
import argparse
import sys
import csv

from tap import Tap
import numpy as np

def f(v: str) -> str:
    """Useful function.

    Docstring is google-style.

    params:
    - v: Interesting input.

    returns: Useful value.
    """
    return v

class Args(Tap):
    cool_value: str = "foo"  # Name of column with cool value.

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

    # Write TSV
    w = csv.DictWriter(sys.stdout, fieldnames=[*reader.fieldnames, "bar"], delimiter="\t")
    w.writeheader()
    for row in reader:
        bar = f(row[args.cool_value])
        w.writerow({**row, "bar": bar})
