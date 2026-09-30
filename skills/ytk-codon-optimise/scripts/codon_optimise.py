"""Rewrite a coding sequence with the most-used codon of every amino acid.

    python3 codon_optimise.py --input genes.csv

The protein is left exactly as it was: translated before and after and compared,
and if it differs at all nothing is kept. The gene keeps its own start and stop
codon. The input file is never touched - an optimised_genes/ folder is written
instead, with the old sequence beside the new one.

**ytk-remove-cut-sites has to run after this.** Rewriting codons can create a
BsmBI or BsaI site, and nothing here looks for one. This exits 1 when a site is
present afterwards, so the next step cannot be forgotten by accident.

Two things get worse, not better, and both are printed rather than fixed:

    GC drops   yeast's preferred codons are AT-rich, so a gene at 50% GC comes
               out near 32%. Low GC makes PCR and synthesis harder.
    repeats    the same peptide motif twice becomes the same DNA twice.

If the gene is going to a synthesis company anyway, their optimiser is better
than this one and free. This is for when the sequence is needed before the order
goes in, or when exactly which codons moved has to be on the record.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import codons       # noqa: E402
import enzymes      # noqa: E402
import metrics      # noqa: E402
import notes        # noqa: E402
import output       # noqa: E402
import sequences    # noqa: E402
import snapgene as sg  # noqa: E402

# Below this, PCR and synthesis both get harder. A guide, not a rule - but the
# rewrite drops GC every time, so the number is worth saying out loud.
GC_FLOOR = 35.0

# More changed codons than this and the table is a wall rather than a record, so
# the console gets the counts and changes.csv gets every row.
PRINT_LIMIT = 40


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", nargs="+", required=True,
                        help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--taxid", type=int, default=codons.YEAST,
                        help=f"the codon usage table to use "
                             f"(default {codons.YEAST}, Saccharomyces cerevisiae)")
    parser.add_argument("--outdir",
                        help="where to write the optimised sequences (default: "
                             "a ytk_output/ folder beside the input)")
    args = parser.parse_args()

    found = sequences.read_all(args.input, args.feature)
    outdir = output.folder(args.outdir, args.input[0])
    # A site cleared by an earlier run cannot be measured off the sequence
    # afterwards, so it travels in its own column and is picked up here.
    was_removed = {}
    for path in args.input:
        was_removed.update(sg.read_column(path, sg.REMOVED_HEADERS))

    method = codons.method(args.taxid)
    folder = outdir / "optimised_genes"
    folder.mkdir(parents=True, exist_ok=True)

    print(f"{'gene':22s} {'bp':>5s} {'codons':>6s} {'CAI':>13s} {'GC%':>13s} "
          f"{'repeat':>13s}  sites after")
    handoff, summary, changed_rows, warnings = [], [], [], []
    for name, sequence, plasmid in found:
        try:
            new, changes = codons.recode(sequence, args.taxid)
        except ValueError as refused:
            sys.exit(f"{name}: {refused}\n\nNothing was changed. Run ytk-cds-qc "
                     f"first and fix what it reports.")

        before, after = metrics.measure(sequence), metrics.measure(new)
        left = enzymes.in_sequence(new)
        # A rewrite can clear a site by accident as easily as it can make one.
        # It counts as removed either way: it cannot be measured off the new
        # sequence, so if it is not recorded here nothing downstream can say the
        # gene ever had it.
        removed = (notes.unpack(was_removed.get(name, ""))
                   | (set(enzymes.in_sequence(sequence)) - set(left)))

        handoff.append([name, plasmid, new, notes.pack(removed), method])
        summary.append([name, sequence, new, "true",
                        notes.summarise(removed, left, method)])
        changed_rows += [[name, change["codon_number"], change["amino_acid"],
                          change["was"], change["now"]]
                         for change in changes]

        print(f"{name[:22]:22s} {len(sequence):5d} {len(changes):6d} "
              f"{before['cai_scer']:5} -> {after['cai_scer']:5} "
              f"{before['gc']:5} -> {after['gc']:5} "
              f"{before['longest_repeat']:5d} -> {after['longest_repeat']:5d}  "
              f"{', '.join(left) if left else 'none'}")

        if after["gc"] < GC_FLOOR:
            warnings.append(f"{name}: GC is now {after['gc']}%, under {GC_FLOOR}%")
        if after["longest_repeat"] > before["longest_repeat"]:
            warnings.append(f"{name}: longest repeat grew from "
                            f"{before['longest_repeat']} to "
                            f"{after['longest_repeat']} bases")
        if changes and len(changes) <= PRINT_LIMIT:
            print("  " + codons.describe(changes).replace("\n", "\n  "))

    output.write_genes_csv(folder / "input_optimised.csv", handoff)
    output.write_summary_csv(folder / "summary.csv", summary)
    with open(folder / "changes.csv", "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "codon_number", "amino_acid", "was", "now"])
        writer.writerows(changed_rows)

    print(f"\n{len(changed_rows)} codons changed across {len(found)} genes. "
          f"Wrote {folder}")
    print(f"Method: {method}. Your input file was not touched.")
    if warnings:
        print("\n" + "\n".join(warnings))
        print("These are not fixed here. Nothing in this plugin fixes them.")

    unsafe = [name for name, _, sequence, *_ in handoff if enzymes.in_sequence(sequence)]
    if unsafe:
        print(f"\nThese now hold a cut site: {', '.join(unsafe)}. Remove them "
              f"before going on:\n"
              f"    python3 \"$CLAUDE_PLUGIN_ROOT/skills/ytk-remove-cut-sites/"
              f"scripts/remove_cut_sites.py\" \\\n"
              f"        --input {folder / 'input_optimised.csv'} --outdir {outdir}")
        sys.exit(1)
    print(f"\nNo cut sites were created. Run ytk-remove-cut-sites anyway if you "
          f"want the check on the record.")


if __name__ == "__main__":
    main()
