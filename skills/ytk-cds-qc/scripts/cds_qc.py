"""Check coding sequences before any cloning is designed.

Two questions: is this a YTK Type 3 coding sequence, and does it hide a BsmBI or
BsaI site that would break Golden Gate?

    python3 cds_qc.py --input genes.csv
    python3 cds_qc.py --sequence ATGAAA...TAA --name my_gene

Nothing is written and nothing is changed. Exits 1 if anything failed.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds          # noqa: E402
import enzymes      # noqa: E402
import report       # noqa: E402
import sequences    # noqa: E402


def collect(inputs, one_sequence, one_name, feature):
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
    parser.add_argument("--sequence", help="one coding sequence")
    parser.add_argument("--name", help="what to call a --sequence")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    args = parser.parse_args()

    if not args.input and not args.sequence:
        parser.error("give --input or --sequence")

    found = collect(args.input, args.sequence, args.name, args.feature)
    log = report.Report()

    print(f"{'sequence':28s} {'bp':>6s} {'codons':>6s} {'BsmBI':>5s} "
          f"{'BsaI':>4s}  verdict")
    passed = 0
    for name, sequence in found:
        faults = cds.problems(sequence)
        notes = cds.warnings(sequence) if not faults else []
        bsmbi = enzymes.count(sequence, "BsmBI")
        bsai = enzymes.count(sequence, "BsaI")
        if faults:
            verdict = "not a Type 3 CDS"
        elif notes:
            verdict = "usable, with warnings"
        else:
            verdict = "ok"
            passed += 1
        codons = len(sequence) // 3 if not len(sequence) % 3 else 0
        print(f"{name[:28]:28s} {len(sequence):6d} "
              f"{codons if codons else '-':>6} {bsmbi:5d} {bsai:4d}  {verdict}")
        for fault in faults:
            log.block(name, fault)
        log.add_warnings(name, notes)

    print(f"\n{passed} of {len(found)} are clean Type 3 coding sequences.")

    table = log.table()
    if table:
        print(f"\n{table}")
    if log.blocked:
        print(f"\n{log.question()}")
        sys.exit(1)
    if log.warned:
        print("\nThese are valid coding sequences, but a cut site inside one has "
              "to go before it will clone: Golden Gate would cut the part apart. "
              "Removing it means changing the DNA, so nothing is changed here.")
        sys.exit(1)


if __name__ == "__main__":
    main()
