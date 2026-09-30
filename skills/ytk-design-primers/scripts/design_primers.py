"""Design PCR primers that amplify a coding sequence with the YTK flanks on.

The product of the PCR is already YTK-compatible, so it drops straight into an
entry vector in a BsmBI Golden Gate reaction.

    python3 design_primers.py --input genes.csv
    python3 design_primers.py --input parts.csv --type 3a

The part type defaults to 3, a whole coding sequence, and is printed on every
run. It is never read off the sequence: 3 against 3a against 3b is a design
decision, not a property of the DNA.

One CSV comes out, ytk_primers.csv, two rows per sequence. Input sequences are
never changed. Without --outdir it goes in a ytk_output/ folder beside the
input.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds             # noqa: E402
import flanks          # noqa: E402
import output          # noqa: E402
import primers         # noqa: E402
import report          # noqa: E402
import sequences       # noqa: E402

# The part type this skill designs for unless --type says otherwise.
DEFAULT_TYPE = "3"

# The spare bases at the outer end of each primer. Four is the shortest that
# still leaves BsmBI room to cut, and nothing has ever wanted a different
# number, so it is a constant rather than a flag.
PAD = flanks.MIN_PAD


def collect(inputs, feature):
    """Sequences from every input file, as (name, sequence)."""
    found = []
    for path in inputs:
        try:
            found += [(name, seq) for name, seq, _ in sequences.read(path, feature)]
        except sequences.AmbiguousCDS as ambiguous:
            sys.exit(f"{ambiguous}\n\n{ambiguous.table()}\n\n"
                     "Say which one with --feature NAME.")
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True,
                        help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--type", default=DEFAULT_TYPE,
                        help=f"YTK part type (default {DEFAULT_TYPE})")
    parser.add_argument("--outdir",
                        help="where to write ytk_primers.csv "
                             "(default: a ytk_output/ folder beside the input)")
    args = parser.parse_args()

    adapter_row = flanks.adapters(args.type)
    found = collect(args.input, args.feature)
    log = report.Report()
    print(f"designing type {args.type} primers\n")

    rows, designs = [], []
    for name, sequence in found:
        faults = cds.problems(sequence, adapter_row["coding"])
        if faults:
            for fault in faults:
                log.block(name, f"does not fit type {args.type}: it {fault}")
            continue
        designed = primers.design(name, sequence, adapter_row, PAD)
        if "blocked" in designed:
            log.block(name, designed["blocked"])
            continue
        log.add_warnings(name, primers.all_warnings(designed))
        designs.append(designed)
        rows += primers.rows_for_csv(designed)

    if designs:
        out = output.folder(args.outdir, args.input[0]) / "ytk_primers.csv"
        with open(out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=primers.CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(designs)} pairs -> {out}\n")
        print(f"{'oligo':32s} {'bp':>3s} {'bind':>4s} {'GC%':>4s} {'Tm':>3s}")
        for row in rows:
            print(f"{row['oligo'][:32]:32s} {row['full_length']:3d} "
                  f"{row['binding_region_length']:4d} {row['bind_region_gc']:4d} "
                  f"{row['tm_phusion']:3d}")
        # An annealing temperature belongs to a pair, not to an oligo, so it is
        # printed here and not in the CSV, where a row is one oligo.
        print(f"\n{'pair':32s} {'anneal at':>9s}")
        for designed in designs:
            print(f"{designed['name'][:32]:32s} {designed['annealing_temp']:8d} C")
        print("\nCheck each pair on the NEB Tm calculator before you order.")

    table = log.table()
    if table:
        print(f"\n{table}")
    question = log.question()
    if question:
        print(f"\n{question}")
        sys.exit(1)


if __name__ == "__main__":
    main()
