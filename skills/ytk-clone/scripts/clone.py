"""Turn flanked YTK fragments into Type 3 part plasmids.

Reads a FASTA or CSV file of fragments, each one already carrying the YTK
flanks, as ordered from a synthesis company. Use ytk-add-overhangs to design
them. For each fragment it writes two SnapGene files: the fragment on its own
(linear) and the finished plasmid (circular).

    python3 clone.py --input fragments.csv --outdir out/

Plasmid names come from the plasmid name column of the CSV, because real
plasmid numbers are not a tidy series. Without that column a plasmid is named
<fragment>_<backbone>, for example Pi_fim_NCS_c1_pYTK001.

Input sequences are never changed. If one looks wrong the script stops and
says so.
"""
import argparse
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import snapgene as sg

sg.require("biopython", "Bio")
sg.require("pydna")

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


def cut_insert(fragment):
    """Cut the fragment with BsmBI and keep the piece between the two designed
    cuts.

    Some genes hold a BsmBI site inside the coding sequence. Cutting then
    gives extra pieces in the middle. Joining them back puts the gene together
    again, so only the two outer cuts count and the gene is left as it was.
    """
    pieces = Dseqrecord(fragment).cut(BsmBI)
    if len(pieces) < 3:
        raise WrongFragment(
            "there is no pair of BsmBI sites here, so this looks like a bare "
            "gene. Run ytk-add-overhangs on it first")
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


def cut_backbone(backbone):
    """Cut the entry vector with BsmBI and keep the larger piece."""
    return max(Dseqrecord(backbone, circular=True).cut(BsmBI), key=len)


def assemble(fragment, backbone):
    """Build the finished circular plasmid.

    A circle has no natural first base, so the map has to pick one. The start
    is put at base 1 of the backbone. That point is always in the backbone, so
    the gene is never split across the start of the file, and the same fragment
    always gives the same map.
    """
    insert = cut_insert(fragment)
    part = str(insert.seq).upper()

    # The cut must take a piece straight out of the fragment. A fragment with a
    # cut site inside it is cut and joined back on the way through, so this is
    # the one place where a change could slip in unseen.
    if part not in fragment.upper():
        raise SequenceChanged(
            "the piece cut out is not in the fragment exactly as it was given")

    junctions = insert_junctions(insert)
    if junctions != part_type_junctions():
        raise WrongFragment(
            f"its overhangs are {junctions[0]} and {junctions[1]}, but a Type "
            f"{PART_TYPE} part for {BACKBONE_NAME} needs "
            f"{part_type_junctions()[0]} and {part_type_junctions()[1]}. "
            "This is the wrong part type for this entry vector")

    plasmid = (cut_backbone(backbone) + insert).looped()
    plasmid = sg.rotate_to(str(plasmid.seq).upper(), backbone[:ANCHOR_LENGTH])

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
    """True when the fragment holds a BsmBI or BsaI site beyond its own flanks.

    A flanked fragment carries two of each by design, one at each end.
    """
    return (len(re.findall("CGTCTC|GAGACG", fragment, re.I)) > 2
            or len(re.findall("GGTCTC|GAGACC", fragment, re.I)) > 2)


# --- main --------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True,
                    help="FASTA or CSV file of flanked fragments")
    ap.add_argument("--outdir", required=True, help="where to write the .dna files")
    args = ap.parse_args()

    genes = sg.read_genes(args.input)
    backbone = backbone_sequence()

    # Build every plasmid first and write nothing yet. If one fragment fails a
    # check, the run stops with no files written, instead of leaving half a
    # batch on disk for someone to find later.
    built = []
    for name, fragment, plasmid_name in genes:
        plasmid_name = plasmid_file_name(name, plasmid_name)
        try:
            built.append((name, fragment, plasmid_name, assemble(fragment, backbone)))
        except WrongFragment as problem:
            sys.exit(f"{name}: {problem}.\n"
                     f"Nothing was written. Your sequence was not changed.")
        except SequenceChanged as problem:
            sys.exit(f"{name}: {problem}\n"
                     f"Nothing was written. Your sequence was not changed.\n"
                     f"Send this fragment to whoever maintains the plugin.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    width = max(len(p) for _, _, p, _ in built) + 2
    print(f"{'fragment':24s} {'plasmid':{width}s} {'frag bp':>8s} "
          f"{'plasmid bp':>11s}  internal site")
    for name, gene, plasmid_name, plasmid in built:
        sg.write_dna(outdir / f"{name}.dna", gene, circular=False,
                     notes_type="Natural")
        sg.write_dna(outdir / f"{plasmid_name}.dna", plasmid, circular=True,
                     notes_type="Synthetic", description="synthetic circular DNA.")

        print(f"{name:24s} {plasmid_name:{width}s} {len(gene):>8d} {len(plasmid):>11d}"
              f"  {'yes' if has_internal_site(gene) else 'no'}")

    flagged = [name for name, gene, _, _ in built if has_internal_site(gene)]
    if flagged:
        print(f"\n{len(flagged)} fragment(s) hold a BsmBI or BsaI site beyond "
              f"their flanks: {', '.join(flagged)}")
        print("They were kept exactly as given. Each one was checked base by "
              "base against the finished plasmid.")
        print("Worth knowing at the bench: those parts cannot be re-cut with "
              "the same enzyme later.")

    print(f"\nwrote {2 * len(genes)} files to {outdir}")


if __name__ == "__main__":
    main()
