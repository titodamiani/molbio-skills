"""Turn gene sequences into YTK Type 3 part plasmids.

Reads a FASTA or CSV file of genes. For each gene it writes two SnapGene
files: the gene on its own (linear) and the finished plasmid (circular).

    python3 clone.py --input genes.csv --outdir out/

Plasmid names come from the plasmid name column of the CSV, because real
plasmid numbers are not a tidy series. Without that column a plasmid is named
<gene>_<backbone>, for example Pi_fim_NCS_c1_pYTK001.

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

from Bio.Restriction import BsmBI          # noqa: E402
from pydna.dseqrecord import Dseqrecord    # noqa: E402

PARTS_TABLE = ROOT / "data" / "ytk_parts.tsv"

# The entry vector every Type 3 part goes into. Its name is also used to name
# plasmid files when the gene file does not give a name.
BACKBONE_NAME = "pYTK001"

# The gene is ordered with these flanks around it. They carry the BsmBI sites
# that cut the part out, and the BsaI sites that put it into a later assembly.
BEGIN_FLANK = "actcgacaacCGTCTCatcGGTCTCaT"
END_FLANK = "ATCCtGAGACCtGAGACGgttgtggtgt"

# How many bases of the backbone are used to fix where the map starts.
ANCHOR_LENGTH = 40

# --- the backbone ------------------------------------------------------


def backbone_sequence():
    """The entry vector sequence, from the shared parts table."""
    with open(PARTS_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["name"] == BACKBONE_NAME:
                return row["sequence"].upper()


def plasmid_file_name(gene_name, given_name):
    """What to call the plasmid file.

    The gene file names it when it can. Otherwise the name says which gene
    went into which vector, so a folder of files still reads clearly later.
    """
    return given_name or f"{gene_name}_{BACKBONE_NAME}"


# --- assembly ----------------------------------------------------------


def cut_insert(gene):
    """Cut the flanked gene with BsmBI and keep the piece between the two
    designed cuts.

    Some genes hold a BsmBI site inside the coding sequence. Cutting then
    gives extra pieces in the middle. Joining them back puts the gene together
    again, so only the two outer cuts count and the gene is left as it was.
    """
    pieces = Dseqrecord(BEGIN_FLANK + gene + END_FLANK).cut(BsmBI)
    kept = pieces[1]
    for piece in pieces[2:-1]:
        kept = kept + piece
    return kept


def cut_backbone(backbone):
    """Cut the entry vector with BsmBI and keep the larger piece."""
    return max(Dseqrecord(backbone, circular=True).cut(BsmBI), key=len)


def assemble(gene, backbone):
    """Build the finished circular plasmid.

    A circle has no natural first base, so the map has to pick one. The start
    is put at base 1 of the backbone. That point is always in the backbone, so
    the gene is never split across the start of the file, and the same gene
    always gives the same map.
    """
    plasmid = (cut_backbone(backbone) + cut_insert(gene)).looped()
    plasmid = sg.rotate_to(str(plasmid.seq).upper(), backbone[:ANCHOR_LENGTH])

    # The gene must come out of the assembly exactly as it went in. Genes with
    # a cut site inside them are cut and joined back on the way through, so
    # this is the one place where a change could slip in unseen. Checking here
    # means a changed sequence stops the run instead of reaching a .dna file.
    if sg.find_in_circle(plasmid, gene) is None:
        raise SequenceChanged(
            "the gene is not in the finished plasmid exactly as it was given")
    return plasmid


class SequenceChanged(Exception):
    """The gene did not survive the assembly unchanged."""


def has_internal_site(gene):
    """True when the gene holds a BsmBI or BsaI site of its own."""
    return bool(re.search("CGTCTC|GAGACG|GGTCTC|GAGACC", gene, re.I))


# --- main --------------------------------------------------------------


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="FASTA or CSV file of genes")
    ap.add_argument("--outdir", required=True, help="where to write the .dna files")
    args = ap.parse_args()

    genes = sg.read_genes(args.input)
    backbone = backbone_sequence()

    # Build every plasmid first and write nothing yet. If one gene fails the
    # check, the run stops with no files written, instead of leaving half a
    # batch on disk for someone to find later.
    built = []
    for name, gene, plasmid_name in genes:
        plasmid_name = plasmid_file_name(name, plasmid_name)
        try:
            built.append((name, gene, plasmid_name, assemble(gene, backbone)))
        except SequenceChanged as problem:
            sys.exit(f"{name}: {problem}\n"
                     f"Nothing was written. Your sequence was not changed.\n"
                     f"Send this gene to whoever maintains the plugin.")

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    width = max(len(p) for _, _, p, _ in built) + 2
    print(f"{'gene':24s} {'plasmid':{width}s} {'gene bp':>8s} "
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
        print(f"\n{len(flagged)} gene(s) hold a BsmBI or BsaI site inside them: "
              f"{', '.join(flagged)}")
        print("They were kept exactly as given. Each one was checked base by "
              "base against the finished plasmid.")
        print("Worth knowing at the bench: those genes cannot be re-cut with "
              "the same enzyme later.")

    print(f"\nwrote {2 * len(genes)} files to {outdir}")


if __name__ == "__main__":
    main()
