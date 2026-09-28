"""Check coding sequences before any cloning is designed.

Two questions: is this a YTK Type 3 coding sequence, and does it hide a BsmBI or
BsaI site that would break Golden Gate?

    python3 cds_qc.py --input genes.csv
    python3 cds_qc.py --sequence ATGAAA...TAA --name my_gene

Nothing is written and nothing is changed. Exits 1 if anything failed.

Asked for by name, and only then, it can also remove an internal cut site by
swapping one codon for another that makes the same amino acid:

    python3 cds_qc.py --input genes.csv --remove-sites --outdir fixed/

That writes a new file and never touches the input. Every changed base is
printed, and the protein is compared before and after.
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
import silent       # noqa: E402


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
    parser.add_argument("--remove-sites", action="store_true",
                        help="swap a codon to remove an internal cut site, "
                             "keeping the protein the same. Writes a new file "
                             "and leaves the input alone")
    parser.add_argument("--outdir", default=".",
                        help="where to write the corrected sequences")
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

    if args.remove_sites:
        return correct(found, log, Path(args.outdir))

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


def correct(found, log, outdir):
    """Remove internal cut sites, and say exactly what changed.

    Only reached when --remove-sites was passed. The input file is never
    touched: a new FASTA is written beside it.
    """
    if log.blocked:
        print(f"\n{log.table()}\n\n{log.question()}")
        print("\nNothing was changed. A sequence has to be a valid Type 3 CDS "
              "before a codon can be swapped safely.")
        sys.exit(1)

    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "corrected.fasta"
    lines, touched, stubborn = [], 0, []
    for name, sequence in found:
        fixed, changes, left = silent.remove_sites(sequence)
        lines.append(f">{name}\n{fixed}")
        if changes:
            touched += 1
            print(f"\n{name}: {len(changes)} base(s) changed")
            print("  " + silent.describe(changes).replace("\n", "\n  "))
        if left:
            stubborn.append(name)
    out.write_text("\n".join(lines) + "\n")

    print(f"\n{touched} of {len(found)} sequences changed. Wrote {out}")
    print("Your input file was not touched.")
    if stubborn:
        print(f"\nStill holding a site, because no synonymous codon removes it: "
              f"{', '.join(stubborn)}. Those need doing by hand.")
        sys.exit(1)


if __name__ == "__main__":
    main()
