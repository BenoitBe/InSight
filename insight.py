#!/usr/bin/env python3
"""InSight — biological reading of a finished Deimos run.

    python insight.py proteomics_output
    python insight.py proteomics_output --out figures --k 6 --volcano A_vs_B --volcano C_vs_D
"""
import argparse
import sys

from insight.run import run_insight


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("deimos_out", help="Deimos output folder (contains ProteomicAnalysis_Results*.xlsx)")
    ap.add_argument("--out", help="output folder (default: <deimos_out>/insight)")
    ap.add_argument("--xlsx", help="results workbook, if not the latest in deimos_out")
    ap.add_argument("--reference-fasta", help="reference proteome of the orthologues (protein names)")
    ap.add_argument("--k", type=int, help="number of response programmes (default: automatic)")
    ap.add_argument("--volcano", action="append", help="contrast to draw (repeatable; default: "
                    "the comparisons with most differential proteins)")
    ap.add_argument("--max-volcanos", type=int, default=4)
    ap.add_argument("--order-regex", help=r"regex reading the ordered level in condition names, "
                    r"e.g. 'T(\d+)' (default: automatic)")
    ap.add_argument("--categorical", action="store_true",
                    help="never read an order in condition names (no trajectories)")
    ap.add_argument("--title", help="report title")
    a = ap.parse_args(argv)
    ov = {}
    if a.order_regex:
        ov["order_regex"] = a.order_regex
    if a.categorical:
        ov["categorical"] = True
    run_insight(a.deimos_out, a.out, a.xlsx, a.reference_fasta, a.k, a.volcano, a.max_volcanos,
                ov or None, a.title)
    return 0


if __name__ == "__main__":
    sys.exit(main())
