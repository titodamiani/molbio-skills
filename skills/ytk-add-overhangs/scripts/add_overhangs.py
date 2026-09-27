"""Put YTK Golden Gate flanks around a sequence, so it can be ordered.

The flanks carry the BsmBI sites that cut the part out of the synthesised
fragment, the BsaI sites that put it into a later assembly, and the four-base
overhang pair that fixes where the part sits in a YTK assembly.

    python3 add_overhangs.py --sequence ATGGCG... --type 3 --name Pi_fim_NCS_c1
    python3 add_overhangs.py --input parts.csv --type 3 --outdir out/

The part type must be given. It is never guessed from the sequence.

One SnapGene file is written per sequence, called <name>_oh.dna. Use
--format fasta or --format csv for plain text instead.

Input sequences are never changed. If one does not fit the type asked for,
the script says so and stops.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import cds        # noqa: E402
import enzymes    # noqa: E402
import flanks     # noqa: E402
import sequences  # noqa: E402
import snapgene as sg

OVERHANG_TABLE = ROOT / "data" / "ytk_overhangs.tsv"

INSERT_COLOR = "#66ccff"


def read_overhangs():
    """The per-type adapter pair, from the shared table."""
    with open(OVERHANG_TABLE, newline="") as fh:
        return {row["part_type"]: row for row in csv.DictReader(fh, delimiter="\t")}


def check(name, sequence, adapters, allow_no_stop):
    """Refuse a sequence that does not fit the part type. Returns warnings.

    A hard stop is for a design that would be wrong DNA: a junction that does
    not form, a reading frame that does not close, a stop codon in the middle
    of a protein fusion. Everything else is reported and left alone.
    """
    if not re.fullmatch(r"[ACGT]+", sequence):
        sys.exit(f"{name}: the sequence holds something other than A, C, G and T")

    coding = adapters["coding"]
    part_type = adapters["part_type"]
    ends_in_stop = sequence[-3:] in cds.STOP_CODONS

    if coding and len(sequence) % 3:
        sys.exit(f"{name}: length {len(sequence)} is not a whole number of codons")
    # A left adapter of a single T relies on the gene's own ATG to complete the
    # TATG overhang, so without that ATG the junction is simply wrong.
    if coding in ("full", "start") and not sequence.startswith("ATG"):
        sys.exit(f"{name}: type {part_type} needs a sequence starting with ATG")
    if coding == "start" and ends_in_stop:
        sys.exit(f"{name}: type {part_type} is the first half of a protein fusion, "
                 "so it must not end with a stop codon")
    if coding in ("full", "end") and not ends_in_stop and not allow_no_stop:
        sys.exit(f"{name}: type {part_type} has no stop codon at the end, and the "
                 f"flank does not add one. Pass --no-stop-codon {name} if that is "
                 "intended")

    return [f"holds a {enzyme} site inside the sequence"
            for enzyme in enzymes.in_sequence(sequence)]


def flank(sequence, adapters):
    """The sequence to order, with both handles at full length."""
    return flanks.flank(sequence, adapters)


def write_snapgene(path, name, sequence, ordered, adapters):
    """Write the flanked sequence as a linear SnapGene file.

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
    sg.write_dna(path, ordered, circular=False, notes_type="Synthetic",
                 description=f"YTK type {adapters['part_type']} part, ready to order.",
                 features=[feature])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequence", help="one DNA sequence")
    parser.add_argument("--input", help=".fa, .fasta, .csv, .gb, .gbk or .dna")
    parser.add_argument("--feature", help="which feature holds the CDS, for map files")
    parser.add_argument("--type", required=True, help="YTK part type, for example 3")
    parser.add_argument("--name", help="name for a single sequence")
    parser.add_argument("--format", default="dna", choices=["dna", "fasta", "csv"],
                       help="output file format (default dna)")
    parser.add_argument("--outdir", default=".", help="where to write the output")
    parser.add_argument("--no-stop-codon", nargs="+", default=[], metavar="NAME",
                       help="names of sequences that end without a stop codon on purpose")
    parser.add_argument("--approve-unsure", action="store_true",
                       help="go ahead with a part type whose overhangs are not confirmed")
    args = parser.parse_args()

    table = read_overhangs()
    if args.type not in table:
        sys.exit(f"unknown part type {args.type}. Known types: {', '.join(table)}")
    adapters = table[args.type]

    if adapters["unsure"] and not args.approve_unsure:
        sys.exit(f"type {args.type}: {adapters['unsure']}\n"
                 "Nothing was written. Pass --approve-unsure to go ahead anyway.")

    if args.input:
        try:
            parts = [(name, seq.upper())
                     for name, seq, _ in sequences.read(Path(args.input), args.feature)]
        except sequences.AmbiguousCDS as ambiguous:
            sys.exit(f"{ambiguous}\n\n{ambiguous.table()}\n\n"
                     "Say which one with --feature NAME.")
    elif args.name:
        parts = [(args.name, args.sequence.upper())]
    else:
        sys.exit("--name is needed, so the output file can be named after the sequence")

    # A typo here would otherwise fail with a message telling the person to
    # pass a flag they did pass.
    unknown = set(args.no_stop_codon) - {name for name, _ in parts}
    if unknown:
        sys.exit(f"--no-stop-codon names a sequence that is not in the input: "
                 f"{', '.join(sorted(unknown))}")

    rows = []
    for name, sequence in parts:
        warnings = check(name, sequence, adapters, name in args.no_stop_codon)
        rows.append((name, sequence, flank(sequence, adapters), warnings))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for name, sequence, ordered, warnings in rows:
        print(f"\n{name}  type {args.type}  {len(sequence)} bp in, {len(ordered)} bp to order")
        print(ordered)
        for warning in warnings:
            print(f"  note: {name} {warning}")

        if args.format == "dna":
            path = outdir / f"{name}_oh.dna"
            write_snapgene(path, name, sequence, ordered, adapters)
        elif args.format == "fasta":
            path = outdir / f"{name}_oh.fasta"
            path.write_text(f">{name}_oh YTK type {args.type}\n{ordered}\n")
        if args.format != "csv":
            print(f"  wrote {path}")

    if args.format == "csv":
        path = outdir / "ytk_order.csv"
        with open(path, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["name", "part_type", "length", "sequence", "notes"])
            for name, _, ordered, warnings in rows:
                # The plain name, not <name>_oh: ytk-clone reads this column
                # and names its plasmid files from it.
                writer.writerow([name, args.type, len(ordered), ordered,
                                 "; ".join(warnings)])
        print(f"\nwrote {path}")

    if adapters["note"]:
        print(f"\ntype {args.type}: {adapters['note']}")


if __name__ == "__main__":
    main()
