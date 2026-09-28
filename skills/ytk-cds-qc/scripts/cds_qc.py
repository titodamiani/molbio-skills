"""Check coding sequences before any cloning is designed.

Two questions: is this a YTK Type 3 coding sequence, and does it hide a BsmBI or
BsaI site that would break Golden Gate?

    python3 cds_qc.py --input genes.csv

Nothing is written and nothing is changed. Exits 1 if anything failed.

Asked for by name, and only then, it can also remove an internal cut site by
swapping one codon for another that makes the same amino acid:

    python3 cds_qc.py --input genes.csv --remove-sites

That writes corrected_genes.csv into the output folder and never touches the
input. Every changed base is printed, and the protein is compared before and
after. Without --outdir the results go in a ytk_output/ folder beside the input.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds          # noqa: E402
import enzymes      # noqa: E402
import notes        # noqa: E402
import output       # noqa: E402
import report       # noqa: E402
import sequences    # noqa: E402
import silent       # noqa: E402


def collect(inputs, feature):
    """Sequences as (name, sequence, plasmid name or None).

    The plasmid name is kept, not dropped: --remove-sites writes it back out,
    so the numbers from the input survive into the maps.
    """
    found = []
    for path in inputs:
        try:
            found += sequences.read(path, feature)
        except sequences.AmbiguousCDS as ambiguous:
            sys.exit(f"{ambiguous}\n\n{ambiguous.table()}\n\n"
                     "Say which one with --feature NAME.")
    return found


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True,
                        help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--remove-sites", action="store_true",
                        help="swap a codon to remove an internal cut site, "
                             "keeping the protein the same. Writes a new file "
                             "and leaves the input alone")
    parser.add_argument("--outdir",
                        help="where to write the corrected sequences "
                             "(default: a ytk_output/ folder beside the input)")
    args = parser.parse_args()

    found = collect(args.input, args.feature)
    log = report.Report()

    print(f"{'sequence':28s} {'bp':>6s} {'codons':>6s} {'BsmBI':>5s} "
          f"{'BsaI':>4s}  verdict")
    passed = 0
    for name, sequence, _ in found:
        faults = cds.problems(sequence)
        # Not `notes`: that is the module this file imports, and shadowing it
        # here would break the first call added to this function.
        remarks = cds.warnings(sequence) if not faults else []
        bsmbi = enzymes.count(sequence, "BsmBI")
        bsai = enzymes.count(sequence, "BsaI")
        if faults:
            verdict = "not a Type 3 CDS"
        elif remarks:
            verdict = "usable, with warnings"
        else:
            verdict = "ok"
            passed += 1
        codons = len(sequence) // 3 if not len(sequence) % 3 else 0
        print(f"{name[:28]:28s} {len(sequence):6d} "
              f"{codons if codons else '-':>6} {bsmbi:5d} {bsai:4d}  {verdict}")
        for fault in faults:
            log.block(name, fault)
        log.add_warnings(name, remarks)

    print(f"\n{passed} of {len(found)} are clean Type 3 coding sequences.")

    if args.remove_sites:
        return correct(found, log, output.folder(args.outdir, args.input[0]))

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
    touched: a new CSV is written beside it.

    CSV and not FASTA, because a FASTA header cannot carry the plasmid name or
    the notes column, and both have to reach the next step. A notes-free
    carrier of corrected sequences is how an order file ends up claiming a gene
    was always clean.

    It sits beside input_genes.csv in the output folder, and its columns are in
    the same order as the fragment summary, so the two read the same way round.
    """
    if log.blocked:
        print(f"\n{log.table()}\n\n{log.question()}")
        print("\nNothing was changed. A sequence has to be a valid Type 3 CDS "
              "before a codon can be swapped safely.")
        sys.exit(1)

    out = outdir / "corrected_genes.csv"
    rows, touched, stubborn = [], 0, []
    for name, sequence, plasmid in found:
        # The enzymes have to be read off before the swap. remove_sites reports
        # which bases moved and which sites are left, never which enzyme it
        # cleared, so the removed set is what was there minus what remains.
        before = {enzyme for _, _, enzyme in silent.sites(sequence)}
        fixed, changes, left = silent.remove_sites(sequence)
        remaining = {enzyme for _, _, enzyme in left}

        # A note for every row, including the untouched ones. A blank note has
        # to mean "nothing known", never "clean".
        rows.append([name, plasmid or "", fixed,
                     notes.after_removal(before, remaining)])

        if changes:
            touched += 1
            print(f"\n{name}: {len(changes)} base(s) changed")
            print("  " + silent.describe(changes).replace("\n", "\n  "))
        if left:
            stubborn.append(name)

    with open(out, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "plasmid", "sequence", "notes"])
        writer.writerows(rows)

    print(f"\n{touched} of {len(found)} sequences changed. Wrote {out}")
    print("Your input file was not touched.")
    if stubborn:
        print(f"\nStill holding a site, because no synonymous codon removes it: "
              f"{', '.join(stubborn)}. Those need doing by hand.")
        sys.exit(1)


if __name__ == "__main__":
    main()
