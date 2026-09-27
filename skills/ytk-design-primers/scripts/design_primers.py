"""Design PCR primers that amplify a coding sequence with the YTK flanks on.

The product of the PCR is already YTK-compatible, so it drops straight into the
pYTK001 entry vector in a BsmBI Golden Gate reaction.

    python3 design_primers.py --input genes.csv --outdir out/
    python3 design_primers.py --sequence ATGAAA...TAA --name my_gene

One CSV comes out, two rows per sequence. Input sequences are never changed.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import flanks          # noqa: E402
import primers         # noqa: E402
import report          # noqa: E402
import sequences       # noqa: E402

OVERHANG_TABLE = ROOT / "data" / "ytk_overhangs.tsv"
PART_TYPE = "3"
STOP_CODONS = ("TAA", "TAG", "TGA")


def adapters():
    """The Type 3 adapter pair, from the shared table."""
    with open(OVERHANG_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["part_type"] == PART_TYPE:
                return row
    sys.exit(f"{OVERHANG_TABLE}: no row for part type {PART_TYPE}")


def not_a_type_3_cds(sequence):
    """Why this is not a full coding sequence, or None when it is fine."""
    if not set(sequence) <= set("ACGT"):
        return "holds something other than A, C, G and T"
    if not sequence.startswith("ATG"):
        return "does not start with ATG"
    if sequence[-3:] not in STOP_CODONS:
        return "does not end with a stop codon"
    if len(sequence) % 3:
        return f"length {len(sequence)} is not a whole number of codons"
    return None


def collect(inputs, one_sequence, one_name, feature):
    """Sequences from wherever they were given, as (name, sequence)."""
    if one_sequence:
        return [(one_name or "sequence", one_sequence.upper())]
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
    parser.add_argument("--input", nargs="*", default=[],
                        help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--sequence", help="one coding sequence, on the command line")
    parser.add_argument("--name", help="what to call a sequence given with --sequence")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--pad", type=int, default=flanks.MIN_PAD,
                        help=f"spare bases at the outer end "
                             f"({flanks.MIN_PAD}-{flanks.FULL_PAD})")
    parser.add_argument("--outdir", default=".")
    args = parser.parse_args()

    if not args.input and not args.sequence:
        parser.error("give --input or --sequence")
    if not flanks.MIN_PAD <= args.pad <= flanks.FULL_PAD:
        parser.error(f"--pad must be {flanks.MIN_PAD} to {flanks.FULL_PAD}; "
                     f"below {flanks.MIN_PAD} leaves BsmBI no room to cut")

    adapter_row = adapters()
    found = collect(args.input, args.sequence, args.name, args.feature)
    log = report.Report()

    rows = []
    for name, sequence in found:
        wrong = not_a_type_3_cds(sequence)
        if wrong:
            log.block(name, f"not a Type 3 CDS: {wrong}")
            continue
        designed = primers.design(name, sequence, adapter_row, args.pad)
        if "blocked" in designed:
            log.block(name, designed["blocked"])
            continue
        log.add_warnings(name, designed["warnings"])
        rows += primers.rows_for_csv(designed)

    if rows:
        outdir = Path(args.outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        out = outdir / "ytk_primers.csv"
        with open(out, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=primers.CSV_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)
        print(f"{len(rows) // 2} pairs -> {out}\n")
        print(f"{'name':32s} {'bp':>3s} {'bind':>4s} {'GC%':>4s} {'Tm':>3s} {'Ta':>3s}")
        for row in rows:
            print(f"{row['name'][:32]:32s} {row['Full length (bp)']:3d} "
                  f"{row['Binding region length (bp)']:4d} {row['GC%']:4d} "
                  f"{row['Tm Phusion (C)']:3d} {row['Combined Tm Phusion (C)']:3d}")

    table = log.table()
    if table:
        print(f"\n{table}")
    question = log.question()
    if question:
        print(f"\n{question}")
        sys.exit(1)


if __name__ == "__main__":
    main()
