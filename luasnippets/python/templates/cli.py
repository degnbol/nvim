"""Title.

Usage details that aren't better described by a flag's help text.
"""
import argparse
import sys
import csv
import numpy as np

def f(v: str) -> str:
    """Useful function.

    params:
    - v: Interesting input.

    returns: Useful value.
    """
    return v

def get_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-v", "--cool-value", default="foo", help="Name of column with cool value.")
    return parser

def main() -> None:
    parser = get_parser()
    args = parser.parse_args()
    
    # Read TSV
    reader = csv.DictReader(sys.stdin, delimiter="\t")
    if reader.fieldnames is None:
        parser.error("No columns found.")
    elif args.cool_value not in reader.fieldnames:
        parser.error(f"Column {args.cool_value} not found among {' '.join(reader.fieldnames)}.")

    # Write TSV
    w = csv.DictWriter(sys.stdout, fieldnames=[*reader.fieldnames, "bar"], delimiter="\t")
    w.writeheader()
    for row in reader:
        bar = f(row[args.cool_value])
        w.writerow({**row, "bar": bar})

