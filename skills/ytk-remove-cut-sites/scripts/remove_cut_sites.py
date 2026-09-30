"""Remove an internal BsmBI or BsaI site from a coding sequence.

A site inside a coding sequence gets cut during Golden Gate, so the part comes
apart. It is removed by swapping one codon for a different codon of the same
amino acid, which leaves the protein exactly as it was.

    python3 remove_cut_sites.py --input genes.csv

**Run this after ytk-codon-optimise, always.** Rewriting a gene's codons can
create a site that was not there before, and the optimiser does not look for
one: removing a site is this script's job and only this script's job, so there
is one implementation of it rather than two.

The input file is never touched. A corrected_genes/ folder is written instead,
with the old sequence beside the new one. Every changed base is printed, and the
protein is translated before and after and compared: if it differs at all,
nothing is kept. Without --outdir the results go in a ytk_output/ folder beside
the input. Exits 1 if a site could not be removed.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds          # noqa: E402
import notes        # noqa: E402
import output       # noqa: E402
import report       # noqa: E402
import sequences    # noqa: E402
import silent       # noqa: E402
import snapgene as sg  # noqa: E402


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True,
                        help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--outdir",
                        help="where to write the corrected sequences (default: "
                             "a ytk_output/ folder beside the input)")
    args = parser.parse_args()

    found = sequences.read_all(args.input, args.feature)
    outdir = output.folder(args.outdir, args.input[0])
    # Whatever an earlier step already did rides on the row with the sequence.
    # Picked up rather than recomputed, because neither fact can be measured off
    # a sequence after the event: a swapped-out site is gone, and the table that
    # rewrote the codons is not written into the DNA.
    was_removed = {}
    optimised = {}
    for path in args.input:
        was_removed.update(sg.read_column(path, sg.REMOVED_HEADERS))
        optimised.update(sg.read_column(path, sg.CODON_OPT_HEADERS))

    log = report.Report()
    for name, sequence, _ in found:
        for fault in cds.problems(sequence, "full"):
            log.block(name, fault)
    if log.blocked:
        print(f"{log.table()}\n\n{log.question()}")
        print("\nNothing was changed. A sequence has to be a valid Type 3 CDS "
              "before a codon can be swapped safely.")
        sys.exit(1)

    folder = outdir / "corrected_genes"
    folder.mkdir(parents=True, exist_ok=True)
    handoff, summary, touched, stubborn = [], [], 0, []
    for name, sequence, plasmid in found:
        # The enzymes have to be read off before the swap. remove_sites reports
        # which bases moved and which sites survive, never which enzyme it
        # cleared, so the removed set is what was there minus what remains.
        before = {enzyme for _, _, enzyme in silent.sites(sequence)}
        fixed, changes, left = silent.remove_sites(sequence)
        remaining = {enzyme for _, _, enzyme in left}
        removed = (before - remaining) | notes.unpack(was_removed.get(name, ""))
        method = optimised.get(name, "")

        handoff.append([name, plasmid, fixed, notes.pack(removed), method])
        summary.append([name, sequence, fixed, str(bool(method)).lower(),
                        notes.summarise(removed, remaining, method or None)])

        if changes:
            touched += 1
            print(f"\n{name}: {len(changes)} base(s) changed")
            print("  " + silent.describe(changes).replace("\n", "\n  "))
        if left:
            stubborn.append(name)

    output.write_genes_csv(folder / "input_corrected.csv", handoff)
    output.write_summary_csv(folder / "summary.csv", summary)

    print(f"\n{touched} of {len(found)} sequences changed. Wrote {folder}")
    print("Your input file was not touched.")
    if stubborn:
        print(f"\nStill holding a site, because no synonymous codon removes it: "
              f"{', '.join(stubborn)}. Those need doing by hand.")
        sys.exit(1)


if __name__ == "__main__":
    main()
