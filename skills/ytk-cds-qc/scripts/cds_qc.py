"""Check coding sequences before any cloning is designed, and measure them.

Two questions: is this a YTK Type 3 coding sequence, and does it hide a BsmBI or
BsaI site that would break Golden Gate?

    python3 cds_qc.py --input genes.csv

Writes input.csv, a copy of what it read, and changes nothing. Whatever format
the input was, that copy is the CSV every later step reads, so a FASTA or a .dna
file does not have to be handled again further down. Exits 1 if anything failed.

**Nothing here can change a sequence.** GC, the longest run of one base, the
longest exact repeat and the yeast CAI are measured and written into input.csv
alongside each sequence, so the numbers say which fixing skill a gene needs:

    a BsmBI or BsaI count above zero  ->  ytk-remove-cut-sites
    a low CAI                         ->  ytk-codon-optimise, perhaps

Without --outdir the results go in a ytk_output/ folder beside the input.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds          # noqa: E402
import metrics      # noqa: E402
import output       # noqa: E402
import report       # noqa: E402
import sequences    # noqa: E402


def collect(inputs, feature):
    """Sequences as (name, sequence, plasmid name or None).

    The plasmid name is kept, not dropped: it goes into input.csv, so the
    numbers from the input survive into the maps.
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
    parser.add_argument("--outdir",
                        help="where to write input.csv (default: a ytk_output/ "
                             "folder beside the input)")
    args = parser.parse_args()

    found = collect(args.input, args.feature)
    outdir = output.folder(args.outdir, args.input[0])
    log = report.Report()

    print(f"{'sequence':24s} {'bp':>5s} {'codons':>6s} {'BsmBI':>5s} {'BsaI':>4s} "
          f"{'GC%':>5s} {'run':>3s} {'rpt':>3s} {'CAI':>5s}  verdict")
    passed = 0
    for name, sequence, _ in found:
        # "full": this skill asks one question, is this a whole coding
        # sequence, so the type is not a choice here.
        faults = cds.problems(sequence, "full")
        # Not `notes`: that is a module name elsewhere in this plugin, and
        # shadowing it here would break the first call added to this function.
        remarks = cds.warnings(sequence, "full") if not faults else []
        measured = metrics.measure(sequence)
        if faults:
            verdict = "not a Type 3 CDS"
        elif remarks:
            verdict = "usable, with warnings"
        else:
            verdict = "ok"
            passed += 1
        codons = len(sequence) // 3 if not len(sequence) % 3 else 0
        cai = measured["cai_scer"]
        print(f"{name[:24]:24s} {len(sequence):5d} "
              f"{codons if codons else '-':>6} {measured['bsmbi']:5d} "
              f"{measured['bsai']:4d} {measured['gc']:5} "
              f"{measured['longest_run']:3d} {measured['longest_repeat']:3d} "
              f"{cai if cai is not None else '-':>5}  {verdict}")
        for fault in faults:
            log.block(name, fault)
        log.add_warnings(name, remarks)

    print(f"\n{passed} of {len(found)} are clean Type 3 coding sequences.")
    print(f"wrote {output.write_input_csv(found, outdir)}")
    print("Your input file was not touched, and nothing here can change a sequence.")

    table = log.table()
    if table:
        print(f"\n{table}")
    if log.blocked:
        print(f"\n{log.question()}")
        sys.exit(1)
    if log.warned:
        print("\nThese are valid coding sequences, but a cut site inside one has "
              "to go before it will clone: Golden Gate would cut the part apart. "
              "Removing it means changing the DNA, so nothing is changed here. "
              "Use ytk-remove-cut-sites, and only if the person asks for it.")
        sys.exit(1)


if __name__ == "__main__":
    main()
