"""Put YTK Golden Gate flanks around a sequence, so it can be ordered.

The flanks carry the BsmBI sites that cut the part out of the synthesised
fragment, the BsaI sites that put it into a later assembly, and the four-base
overhang pair that fixes where the part sits in a YTK assembly.

    python3 add_overhangs.py --input parts.csv
    python3 add_overhangs.py --input parts.csv --type 3a --outdir out/

The part type defaults to 3, a whole coding sequence, and is printed on every
run. Types 1 to 8b all work. It is still never read off the sequence: 3 against 3a
against 3b is a design decision, not a property of the DNA.

Everything lands in a fragments/ folder under the output folder: one GenBank map
per sequence, called <name>.gb, plus summary.csv. That summary is both the file a
synthesis order is placed from and the input to ytk-clone. Use --format dna for
SnapGene files instead. Without --outdir the results go in a ytk_output/ folder
beside the input.

Input sequences are never changed. If one does not fit the type asked for,
the script says so and stops.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds        # noqa: E402
import enzymes    # noqa: E402
import flanks     # noqa: E402
import notes      # noqa: E402
import output     # noqa: E402
import sequences  # noqa: E402
import snapgene as sg

# A whole coding sequence. It is a default and not a guess - it is printed on
# every run, and --type overrides it.
PART_TYPE = "3"

INSERT_COLOR = "#66ccff"


def check(name, sequence, adapters, allow_no_stop):
    """Refuse a sequence that does not fit the part type.

    A hard stop is for a design that would be wrong DNA: a junction that does
    not form, a reading frame that does not close, a stop codon in the middle
    of a protein fusion. The checks themselves live in lib/cds.py and are the
    same ones ytk-design-primers uses.

    Every fault is reported at once, so one run tells you everything to fix.

    A cut site inside the sequence is not checked here. It is measured once,
    for the notes column, and printed from there, so the screen and the order
    file cannot disagree about it.
    """
    faults = cds.problems(sequence, adapters["coding"])
    if allow_no_stop:
        faults = [fault for fault in faults if fault != cds.NO_STOP_CODON]
    if not faults:
        return
    message = "\n".join(f"{name}: type {adapters['part_type']} {fault}"
                        for fault in faults)
    if cds.NO_STOP_CODON in faults:
        message += (f"\nThe flank does not add one. Pass --no-stop-codon {name} "
                    f"if that is intended.")
    sys.exit(message)


def write_labelled_map(path, name, sequence, ordered, adapters):
    """Write the flanked sequence as a linear map, SnapGene .dna or GenBank.

    Only the sequence itself is labelled. The cut sites are left off on
    purpose: SnapGene shows those live under Enzymes, so a fixed label there
    would only go stale.
    """
    start = flanks.insert_offset(adapters)
    feature = {
        "name": name,
        "type": "CDS" if adapters["coding"] else "misc_feature",
        "start": start,
        "end": start + len(sequence),
        "strand": 1,
        "color": INSERT_COLOR,
        "wrap_end": 0,
    }
    sg.write_map(path, ordered, circular=False, notes_type="Synthetic",
                 description=f"YTK type {adapters['part_type']} part, ready to order",
                 features=[feature])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True,
                        help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--type", default=PART_TYPE,
                       help=f"YTK part type (default {PART_TYPE})")
    parser.add_argument("--format", default="genbank", choices=["genbank", "dna"],
                       help="format of the per-fragment maps (default genbank, "
                            "which SnapGene also opens)")
    parser.add_argument("--outdir",
                        help="where to write the output "
                             "(default: a ytk_output/ folder beside the input)")
    parser.add_argument("--no-stop-codon", nargs="+", default=[], metavar="NAME",
                       help="names of sequences that end without a stop codon on purpose")
    args = parser.parse_args()

    adapters = flanks.adapters(args.type)

    try:
        parts = [(name, seq.upper(), plasmid) for name, seq, plasmid
                 in sequences.read(Path(args.input), args.feature)]
    except sequences.AmbiguousCDS as ambiguous:
        sys.exit(f"{ambiguous}\n\n{ambiguous.table()}\n\n"
                 "Say which one with --feature NAME.")

    # A typo here would otherwise fail with a message telling the person to
    # pass a flag they did pass.
    unknown = set(args.no_stop_codon) - {name for name, _, _ in parts}
    if unknown:
        sys.exit(f"--no-stop-codon names a sequence that is not in the input: "
                 f"{', '.join(sorted(unknown))}")

    # What cds_qc did to each gene, read as facts and not as prose. Which
    # enzymes it cleared cannot be measured after the swap, so that one is
    # carried; what is still in the gene is measured again below. Both ride on
    # the same row as the sequence they describe, so neither can end up against
    # the wrong gene.
    removed = sg.read_column(args.input, sg.REMOVED_HEADERS)
    optimised = sg.read_column(args.input, sg.CODON_OPT_HEADERS)

    rows = []
    for name, sequence, _ in parts:
        check(name, sequence, adapters, name in args.no_stop_codon)
        method = optimised.get(name, "")
        note = notes.summarise(notes.unpack(removed.get(name, "")),
                               enzymes.in_sequence(sequence), method)
        rows.append((name, sequence, flanks.flank(sequence, adapters),
                     bool(method), note))

    fragments = output.folder(args.outdir, args.input) / "fragments"
    fragments.mkdir(parents=True, exist_ok=True)
    suffix = ".gb" if args.format == "genbank" else ".dna"

    for name, sequence, ordered, _, note in rows:
        print(f"\n{name}  type {args.type}  {len(sequence)} bp in, {len(ordered)} bp to order")
        print(ordered)
        # The note as it will appear in the order file, word for word, so what
        # is on screen and what gets ordered cannot drift apart. Nothing was
        # done and nothing is left, so there is nothing to print.
        if note:
            print(f"  note: {note}")
        path = fragments / f"{name}{suffix}"
        write_labelled_map(path, name, sequence, ordered, adapters)
        print(f"  wrote {path}")

    # summary.csv is both the synthesis order and the handoff to ytk-clone. The
    # sequence column is the flanked sequence, which is what you paste into an
    # order form. The name column is the plain gene name, not <name>_oh:
    # ytk-clone names its output files from it. Plasmid names are not here -
    # they belong with the input, and ytk-clone reads them from there.
    summary = fragments / "summary.csv"
    with open(summary, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "sequence", "part_type", "codon_opt", "notes"])
        for name, _, ordered, codon_opt, note in rows:
            writer.writerow([name, ordered, args.type,
                             str(codon_opt).lower(), note])
    print(f"\nwrote {summary}")

    if adapters["note"]:
        print(f"\ntype {args.type}: {adapters['note']}")


if __name__ == "__main__":
    main()
