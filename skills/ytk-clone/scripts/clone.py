"""Simulate cloning a flanked fragment into a backbone, and write the map.

This is the in-silico version of what SnapGene's cloning simulation does: cut
the fragment and the backbone with a Type IIS enzyme, check the sticky ends
match, join them, and write the finished map.

Reads fragments that already carry their flanks. Use ytk-add-overhangs to design
them. One circular map per plasmid goes into a <backbone>_maps/ folder, with a
summary.csv beside them. --format picks SnapGene .dna or GenBank.

    python3 clone.py --input fragments/summary.csv
    python3 clone.py --input fragments/summary.csv --outdir out/ \
        --backbone my_vector.gb --enzyme BsaI

Without --outdir the maps go in a ytk_output/ folder beside the input.

The backbone defaults to pYTK001 and the enzyme to BsmBI, which is the YTK entry
reaction. Any Type IIS enzyme Biopython knows works.

Plasmid names come from the plasmid name column of the CSV, because real
plasmid numbers are not a tidy series. Without that column a plasmid is named
<fragment>_<backbone>, for example Pi_fim_NCS_c1_pYTK001.

Input sequences are never changed. If one looks wrong the script stops and
says so.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import snapgene as sg

import deps                               # noqa: E402
import enzymes                           # noqa: E402
import output                             # noqa: E402
import sequences                          # noqa: E402

deps.require("Bio", "pydna")

from Bio import Restriction                # noqa: E402
from Bio.Restriction import BsaI, BsmBI    # noqa: E402
from Bio.Seq import Seq                    # noqa: E402
from pydna.dseqrecord import Dseqrecord    # noqa: E402

PARTS_TABLE = ROOT / "data" / "ytk_parts.tsv"

# The entry vector every Type 3 part goes into. Its name is also used to name
# plasmid files when the fragment file does not give a name.
BACKBONE_NAME = "pYTK001"

# The only part type that belongs in this entry vector.
PART_TYPE = "3"

# How many bases of the backbone are used to fix where the map starts.
ANCHOR_LENGTH = 40

# --- the backbone ------------------------------------------------------


def named_enzyme(name):
    """A Type IIS enzyme by name, from Biopython's set."""
    enzyme = getattr(Restriction, name, None)
    if enzyme is None:
        sys.exit(f"unknown enzyme {name}")
    return enzyme


def load_backbone(path):
    """A backbone from a file, or pYTK001 when none is given."""
    if path is None:
        return BACKBONE_NAME, backbone_sequence()
    try:
        found = sequences.read(path)
    except sequences.AmbiguousCDS:
        # A backbone is a whole vector, so its features do not matter here.
        found = None
    if found:
        return Path(path).stem, found[0][1]
    # Only a map file gets here, because sequences.read raises AmbiguousCDS
    # when a map holds more than one candidate feature.
    return Path(path).stem, sg.read_map(path)["sequence"]


def backbone_sequence():
    """The entry vector sequence, from the shared parts table."""
    with open(PARTS_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["name"] == BACKBONE_NAME:
                return row["sequence"].upper()


def part_type_junctions():
    """The overhang pair a Type 3 part must have, from the shared parts table."""
    with open(PARTS_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["part_type"] == PART_TYPE:
                return row["junction_5"].upper(), row["junction_3"].upper()


def plasmid_file_name(fragment_name, given_name):
    """What to call the plasmid file.

    The fragment file names it when it can. Otherwise the name says which part
    went into which vector, so a folder of files still reads clearly later.
    """
    return given_name or f"{fragment_name}_{BACKBONE_NAME}"





# --- assembly ----------------------------------------------------------


def cut_insert(fragment, enzyme=BsmBI):
    """Cut the fragment with the entry enzyme and keep the piece between the two
    designed cuts.

    Some genes hold a BsmBI site inside the coding sequence. Cutting then
    gives extra pieces in the middle. Joining them back puts the gene together
    again, so only the two outer cuts count and the gene is left as it was.
    """
    pieces = Dseqrecord(fragment).cut(enzyme)
    if len(pieces) < 3:
        raise WrongFragment(
            f"there is no pair of {enzyme} sites here, so this looks like a bare "
            f"gene. Run ytk-add-overhangs on it first")
    kept = pieces[1]
    for piece in pieces[2:-1]:
        kept = kept + piece
    return kept


def insert_junctions(insert):
    """The overhang pair the insert would join a YTK assembly by.

    These come from the inner BsaI sites, not the outer BsmBI ones. The BsmBI
    cut gives the same pair of ends on every part type, because that pair is
    what fits the entry vector. It is the BsaI pair that says which slot of an
    assembly the part belongs in.
    """
    pieces = Dseqrecord(str(insert.seq)).cut(BsaI)
    if len(pieces) < 3:
        raise WrongFragment("there is no pair of BsaI sites inside the fragment")
    five = pieces[1].seq.five_prime_end()[1].upper()
    three = str(Seq(pieces[-2].seq.three_prime_end()[1]).reverse_complement()).upper()
    return five, three


def cut_backbone(backbone, enzyme=BsmBI):
    """Cut the backbone and keep the larger piece.

    In an entry vector the smaller piece is the dropout that the part replaces.
    """
    pieces = Dseqrecord(backbone, circular=True).cut(enzyme)
    if len(pieces) < 2:
        raise WrongFragment(
            f"the backbone has no pair of {enzyme} sites, so there is nothing "
            f"to drop out")
    return max(pieces, key=len)


def assemble(fragment, backbone, enzyme=BsmBI, expect_junctions=True):
    """Build the finished circular plasmid.

    A circle has no natural first base, so the map has to pick one. The start
    goes at base 1 of the backbone, which is always in the backbone, so the gene
    is never split across the start of the file and the same fragment always
    gives the same map.

    Base 1 of a backbone read straight from the kit can sit inside the piece
    that drops out, and then it is gone from the finished plasmid. In that case
    the map starts at the first base of the piece that was kept, which exists by
    definition.
    """
    insert = cut_insert(fragment, enzyme)
    part = str(insert.seq).upper()

    # The cut must take a piece straight out of the fragment. A fragment with a
    # cut site inside it is cut and joined back on the way through, so this is
    # the one place where a change could slip in unseen.
    if part not in fragment.upper():
        raise SequenceChanged(
            "the piece cut out is not in the fragment exactly as it was given")

    junctions = insert_junctions(insert) if expect_junctions else None
    if expect_junctions and junctions != part_type_junctions():
        raise WrongFragment(
            f"its overhangs are {junctions[0]} and {junctions[1]}, but a Type "
            f"{PART_TYPE} part for {BACKBONE_NAME} needs "
            f"{part_type_junctions()[0]} and {part_type_junctions()[1]}. "
            "This is the wrong part type for this entry vector")

    kept = cut_backbone(backbone, enzyme)
    plasmid = (kept + insert).looped()
    anchor = backbone[:ANCHOR_LENGTH].upper()
    circle = str(plasmid.seq).upper()
    if sg.find_in_circle(circle, anchor) is None:
        anchor = str(kept.seq).upper()[:ANCHOR_LENGTH]
    plasmid = sg.rotate_to(circle, anchor)

    # And it must still be there, whole, once the loop is closed and turned.
    # Checking here means a changed sequence stops the run instead of reaching
    # a .dna file.
    if sg.find_in_circle(plasmid, part) is None:
        raise SequenceChanged(
            "the part is not in the finished plasmid exactly as it was cut")
    return plasmid


class SequenceChanged(Exception):
    """The sequence did not survive the assembly unchanged."""


class WrongFragment(Exception):
    """The fragment is not a flanked Type 3 part for this entry vector."""


def has_internal_site(fragment):
    """True when the fragment holds a cut site beyond its own flanks.

    A flanked fragment carries two sites per enzyme by design, one at each end.
    """
    return bool(enzymes.beyond_flanks(fragment))


# --- main --------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", nargs="+", required=True,
                    help="files of flanked fragments (.fa, .fasta, .csv, .gb, "
                         ".gbk, .dna)")
    ap.add_argument("--feature", help="which feature holds the fragment, for map files")
    ap.add_argument("--outdir",
                    help="where to write the maps "
                         "(default: a ytk_output/ folder beside the input)")
    ap.add_argument("--backbone", default=None,
                    help=f"vector to clone into (.dna, .gb, .gbk, FASTA); "
                         f"default {BACKBONE_NAME}")
    ap.add_argument("--enzyme", default="BsmBI",
                    help="Type IIS enzyme to cut with (default BsmBI)")
    ap.add_argument("--format", default="genbank", choices=["genbank", "dna"],
                    help="output map format (default genbank, which SnapGene "
                         "also opens)")
    args = ap.parse_args()

    genes = []
    for path in args.input:
        try:
            genes += sequences.read(path, args.feature)
        except sequences.AmbiguousCDS as ambiguous:
            sys.exit(f"{ambiguous}\n\n{ambiguous.table()}\n\n"
                     "Say which one with --feature NAME.")
    enzyme = named_enzyme(args.enzyme)
    backbone_name, backbone = load_backbone(args.backbone)
    # The Type 3 junction check only means something for the YTK entry reaction.
    standard = args.backbone is None and args.enzyme == "BsmBI"
    if not standard:
        print(f"cloning into {backbone_name} with {args.enzyme}\n")

    # Build every plasmid first and write nothing yet. If one fragment fails a
    # check, the run stops with no files written, instead of leaving half a
    # batch on disk for someone to find later.
    built = []
    for name, fragment, plasmid_name in genes:
        plasmid_name = plasmid_file_name(name, plasmid_name)
        try:
            built.append((name, fragment, plasmid_name,
                          assemble(fragment, backbone, enzyme, standard)))
        except WrongFragment as problem:
            sys.exit(f"{name}: {problem}.\n"
                     f"Nothing was written. Your sequence was not changed.")
        except SequenceChanged as problem:
            sys.exit(f"{name}: {problem}\n"
                     f"Nothing was written. Your sequence was not changed.\n"
                     f"Send this fragment to whoever maintains the plugin.")

    # One folder, named after the backbone actually used, so a run against a
    # different vector cannot be mistaken for an entry-vector run.
    maps = output.folder(args.outdir, args.input[0]) / f"{backbone_name}_maps"
    maps.mkdir(parents=True, exist_ok=True)

    # part_type and notes ride along on the fragment file when ytk-add-overhangs
    # wrote it. Both are optional: a bare FASTA of fragments still clones.
    # Every input file is read, not just the first: with two fragment files the
    # second one's genes would otherwise lose the record of what was cleared
    # out of them, and a blank note is meant to mean "nothing known".
    part_types, notes = {}, {}
    for path in args.input:
        part_types.update(sg.read_column(path, sg.PART_TYPE_HEADERS))
        notes.update(sg.read_column(path, sg.NOTES_HEADERS))

    suffix = ".gb" if args.format == "genbank" else ".dna"
    width = max(len(p) for _, _, p, _ in built) + 2
    print(f"{'fragment':24s} {'plasmid':{width}s} {'frag bp':>8s} "
          f"{'plasmid bp':>11s}  internal site")
    summary = [["name", "part_type", "plasmid", "notes"]]
    for name, gene, plasmid_name, plasmid in built:
        sg.write_map(maps / f"{plasmid_name}{suffix}", plasmid, circular=True,
                     notes_type="Synthetic", description="synthetic circular DNA")
        summary.append([name, part_types.get(name, ""), plasmid_name,
                        notes.get(name, "")])

        print(f"{name:24s} {plasmid_name:{width}s} {len(gene):>8d} {len(plasmid):>11d}"
              f"  {'yes' if has_internal_site(gene) else 'no'}")

    with open(maps / "summary.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(summary)

    # The column above is the report. ytk-cds-qc names these genes first and the
    # notes column carries them into the order file, so all that is left to add
    # is the one thing neither of those says: what it means at the bench.
    if any(has_internal_site(gene) for _, gene, _, _ in built):
        print("\nA fragment marked yes cannot be re-cut with the same enzyme "
              "later. It was kept exactly as given.")

    print(f"\nwrote {len(genes)} maps and summary.csv to {maps}")


if __name__ == "__main__":
    main()
