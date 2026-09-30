"""Simulate cloning a flanked fragment into a backbone, and write the map.

This is the in-silico version of what SnapGene's cloning simulation does: cut
the fragment and the backbone with a Type IIS enzyme, check the sticky ends
match, join them, and write the finished map.

Reads fragments that already carry their flanks. Use ytk-add-overhangs to design
them. One circular map per plasmid goes into plasmids/<backbone>/, with a
summary.csv beside them.

    python3 clone.py --input fragments/summary.csv --genes input.csv
    python3 clone.py --input fragments/summary.csv --outdir out/ \
        --type 5 --backbone pYTK047 --enzyme BsaI

Without --outdir the maps go in a ytk_output/ folder beside the input.

The backbone defaults to pYTK001 and the enzyme to BsmBI, which is the YTK entry
reaction. --backbone takes a published pYTKnnn name or a file of your own. Any
Type IIS enzyme Biopython knows works.

--type says which part type the fragments are. It is a declaration, never read
off the DNA, and it defaults to 3. Every run checks the fragment's overhangs
against the published pair for that type, so a wrong --type stops the run. That
check asks whether the fragment is the type you declared. Whether the part fits
the backbone is a separate question, answered when the two pieces are joined.

The map format is not asked for here. It is read off the fragment maps beside
the input, so the maps come out in whatever format the fragments were written
in. --format overrides that.

Plasmid names come from the plasmid column of --genes, because real plasmid
numbers are not a tidy series. Without one a plasmid is named
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
import flanks                             # noqa: E402
import notes                              # noqa: E402
import output                             # noqa: E402
import sequences                          # noqa: E402

deps.require("Bio", "pydna")

from Bio import Restriction                # noqa: E402
from Bio.Restriction import BsaI, BsmBI    # noqa: E402
from Bio.Seq import Seq                    # noqa: E402
from pydna.dseqrecord import Dseqrecord    # noqa: E402

PARTS_TABLE = ROOT / "data" / "ytk_parts.tsv"

# The 96 published YTK plasmid maps, so --backbone takes a bare pYTKnnn.
PLASMIDS = ROOT / "reference" / "ytk_plasmids"

# The entry vector, used when --backbone says nothing.
BACKBONE_NAME = "pYTK001"

# The part type used when --type says nothing.
DEFAULT_TYPE = "3"

# The flanks always put BsaI inside BsmBI, whatever part type is being built,
# so the inner cut is a fact about the fragment and not a flag. See lib/flanks.
INNER_ENZYME = BsaI

# How many bases of the backbone are used to fix where the map starts.
ANCHOR_LENGTH = 40

# --- the backbone ------------------------------------------------------


def named_enzyme(name):
    """A Type IIS enzyme by name, from Biopython's set."""
    enzyme = getattr(Restriction, name, None)
    if enzyme is None:
        sys.exit(f"unknown enzyme {name}")
    return enzyme


def whole_sequence(path):
    """The whole sequence in a file.

    A backbone is a whole vector, so it must never go through the feature
    picker in lib/sequences: that would hand back one gene out of the vector
    and call it the backbone.
    """
    if Path(path).suffix.lower() in (".fa", ".fasta"):
        return sg.read_genes(path)[0][1]
    return sg.read_map(path)["sequence"]


def load_backbone(given):
    """A backbone by file path or by published name, or pYTK001 by default.

    A path that exists always wins. Otherwise the name is looked for among the
    published maps, so --backbone pYTK047 works with no path typed.
    """
    if given is None or given == BACKBONE_NAME:
        # pYTK001 comes from data/ytk_parts.tsv, stored already turned to
        # SnapGene's own start point. See NOTES.md.
        return BACKBONE_NAME, backbone_sequence()
    path = Path(given)
    if not path.exists():
        path = PLASMIDS / f"{given}.gb"
    if not path.exists():
        sys.exit(f"--backbone {given}: no such file, and no {given}.gb in {PLASMIDS}")
    return path.stem, whole_sequence(path)


def backbone_sequence():
    """The entry vector sequence, from the shared parts table."""
    with open(PARTS_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["name"] == BACKBONE_NAME:
                return row["sequence"].upper()


def plasmid_file_name(fragment_name, given_name, backbone_name):
    """What to call the plasmid file.

    The input file names it when it can. Otherwise the name says which part
    went into which vector, so a folder of files still reads clearly later.
    """
    return given_name or f"{fragment_name}_{backbone_name}"


def map_format(fragment_file, chosen):
    """Which format to write the maps in.

    The format is decided once, when the fragments are written, so reading it
    back off those files means the maps cannot come out in a different format
    from the fragments they were built from. Nobody has to pass the same flag to
    two scripts and remember to keep them the same.

    GenBank when there are no fragment maps to look at, which is the case for a
    bare FASTA of fragments.
    """
    if chosen:
        return chosen
    beside = Path(fragment_file).parent
    return "dna" if any(beside.glob("*.dna")) else "genbank"



# --- assembly ----------------------------------------------------------


def cut_insert(fragment, enzyme):
    """Cut the fragment and keep the piece between the two designed cuts.

    Some genes hold a site for that enzyme inside the coding sequence. Cutting then
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


def insert_junctions(insert, enzyme):
    """The overhang pair the insert would join a YTK assembly by.

    These come from the inner BsaI sites, not the outer BsmBI ones. The outer
    cut gives the same pair of ends on every part type, because that pair is
    what fits the entry vector. It is the BsaI pair that says which slot of an
    assembly the part belongs in.

    A Type IIS site sits outside its own cut, so cutting with BsaI takes both
    BsaI sites away with the ends. Then the piece in hand already has the
    part-type overhangs and there is nothing left to cut.
    """
    if enzyme is INNER_ENZYME:
        pieces = [insert]
    else:
        pieces = Dseqrecord(str(insert.seq)).cut(INNER_ENZYME)
        if len(pieces) < 3:
            raise WrongFragment("there is no pair of BsaI sites inside the fragment")
        pieces = pieces[1:-1]
    five = pieces[0].seq.five_prime_end()[1].upper()
    three = str(Seq(pieces[-1].seq.three_prime_end()[1]).reverse_complement()).upper()
    return five, three


def cut_backbone(backbone, enzyme):
    """Cut the backbone and keep the larger piece.

    In an entry vector the smaller piece is the dropout that the part replaces.
    """
    pieces = Dseqrecord(backbone, circular=True).cut(enzyme)
    if len(pieces) < 2:
        raise WrongFragment(
            f"the backbone has no pair of {enzyme} sites, so there is nothing "
            f"to drop out")
    return max(pieces, key=len)


def accepted_junctions(kept):
    """The overhang pair a cut backbone will take a part by.

    The backbone's 3' end meets the part's 5' overhang, so the pair is read in
    the part's orientation. That way it can be compared with, and printed
    beside, what insert_junctions gives.
    """
    five = kept.seq.five_prime_end()[1].upper()
    three = str(Seq(kept.seq.three_prime_end()[1]).reverse_complement()).upper()
    return three, five


def assemble(fragment, backbone, enzyme, part_type):
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

    # This asks one question: is the fragment the part type that was declared.
    # Whether the part fits this backbone is a different question, and the
    # ligation below answers it.
    junctions = insert_junctions(insert, enzyme)
    expected = flanks.junctions(part_type)
    if junctions != expected:
        raise WrongFragment(
            f"its overhangs are {junctions[0]} and {junctions[1]}, but a type "
            f"{part_type} part needs {expected[0]} and {expected[1]}")

    kept = cut_backbone(backbone, enzyme)
    try:
        plasmid = (kept + insert).looped()
    except TypeError:
        accepts = accepted_junctions(kept)
        raise WrongFragment(
            f"the part offers {junctions[0]} and {junctions[1]}, but the "
            f"backbone accepts {accepts[0]} and {accepts[1]}")
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
    """The fragment is not a flanked part of the type declared, or it does not
    fit the backbone."""


# --- main --------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", nargs="+", required=True,
                    help="files of flanked fragments (.fa, .fasta, .csv, .gb, "
                         ".gbk, .dna)")
    ap.add_argument("--feature", help="which feature holds the fragment, for map files")
    ap.add_argument("--genes",
                    help="the input CSV, for the plasmid names and for what was "
                         "done to each gene")
    ap.add_argument("--outdir",
                    help="where to write the maps "
                         "(default: a ytk_output/ folder beside the input)")
    ap.add_argument("--type", default=DEFAULT_TYPE,
                    help=f"YTK part type the fragments are (default {DEFAULT_TYPE})")
    ap.add_argument("--backbone", default=None,
                    help=f"vector to clone into: a published pYTKnnn name, or a "
                         f"file (.dna, .gb, .gbk, FASTA); default {BACKBONE_NAME}")
    ap.add_argument("--enzyme", default="BsmBI",
                    help="Type IIS enzyme to cut with (default BsmBI)")
    ap.add_argument("--format", default=None, choices=["genbank", "dna"],
                    help="output map format (default: match the fragment maps "
                         "beside the input, or genbank if there are none)")
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
    print(f"cloning type {args.type} fragments into {backbone_name} "
          f"with {args.enzyme}\n")

    # part_type rides along on the fragment file when ytk-add-overhangs wrote
    # it, and is optional: a bare FASTA of fragments still clones. Every input
    # file is read, not just the first, so a second fragment file's genes keep
    # their part type too. A row that disagrees with --type stops the run here,
    # before any plasmid is built.
    part_types = {}
    for path in args.input:
        part_types.update(sg.read_column(path, sg.PART_TYPE_HEADERS))
    for name, written in part_types.items():
        if written and written != args.type:
            sys.exit(f"{name}: the fragment file says part type {written}, but "
                     f"--type says {args.type}. Nothing was written.")

    # Plasmid names, and what was done to each gene, come from the input file.
    # They are facts about the gene and not about the fragment, so they live
    # with the input rather than being copied along the chain.
    def from_genes(headers):
        return sg.read_column(args.genes, headers) if args.genes else {}

    given_names = from_genes(sg.PLASMID_HEADERS)
    removed = from_genes(sg.REMOVED_HEADERS)
    optimised = from_genes(sg.CODON_OPT_HEADERS)

    # Build every plasmid first and write nothing yet. If one fragment fails a
    # check, the run stops with no files written, instead of leaving half a
    # batch on disk for someone to find later.
    built = []
    for name, fragment, plasmid_name in genes:
        plasmid_name = plasmid_file_name(name, given_names.get(name) or plasmid_name,
                                         backbone_name)
        try:
            built.append((name, fragment, plasmid_name,
                          assemble(fragment, backbone, enzyme, args.type)))
        except WrongFragment as problem:
            sys.exit(f"{name}: {problem}.\n"
                     f"Nothing was written. Your sequence was not changed.")
        except SequenceChanged as problem:
            sys.exit(f"{name}: {problem}\n"
                     f"Nothing was written. Your sequence was not changed.\n"
                     f"Send this fragment to whoever maintains the plugin.")

    # One folder per backbone under plasmids/, so a run against a different
    # vector cannot be mistaken for an entry-vector run, and so the tree does
    # not have to change when a second backbone is supported.
    maps = output.folder(args.outdir, args.input[0]) / "plasmids" / backbone_name
    maps.mkdir(parents=True, exist_ok=True)

    suffix = ".gb" if map_format(args.input[0], args.format) == "genbank" else ".dna"
    width = max(len(p) for _, _, p, _ in built) + 2
    print(f"{'fragment':24s} {'plasmid':{width}s} {'frag bp':>8s} "
          f"{'plasmid bp':>11s}  internal BsmBI/BsaI")
    summary = [["plasmid", "sequence", "part_type", "codon_opt", "notes"]]
    flagged = []
    for name, gene, plasmid_name, plasmid in built:
        sg.write_map(maps / f"{plasmid_name}{suffix}", plasmid, circular=True,
                     notes_type="Synthetic", description="synthetic circular DNA")
        # What is still in the gene is measured here, on the fragment in hand,
        # rather than carried from the last step. A flanked fragment has two
        # sites per enzyme by design, so only the extras count.
        extra = enzymes.beyond_flanks(gene)
        if extra:
            flagged.append(name)
        method = optimised.get(name, "")
        summary.append([plasmid_name, name, part_types.get(name, ""),
                        str(bool(method)).lower(),
                        notes.summarise(notes.unpack(removed.get(name, "")),
                                        extra, method)])

        print(f"{name:24s} {plasmid_name:{width}s} {len(gene):>8d} {len(plasmid):>11d}"
              f"  {'yes' if extra else 'no'}")

    with open(maps / "summary.csv", "w", newline="") as fh:
        csv.writer(fh).writerows(summary)

    # The column above is the report. ytk-cds-qc names these genes first and the
    # notes column carries them into the order file, so all that is left to add
    # is the one thing neither of those says: what it means at the bench.
    if flagged:
        print("\nA fragment marked yes cannot be re-cut with the same enzyme "
              "later. It was kept exactly as given.")

    print(f"\nwrote {len(genes)} maps and summary.csv to {maps}")


if __name__ == "__main__":
    main()
