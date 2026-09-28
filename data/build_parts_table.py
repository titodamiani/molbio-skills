"""Build data/ytk_parts.tsv from the published MoClo Yeast Toolkit files.

Source (supplementary material of):
    Lee ME, DeLoache WC, Cervantes B, Dueber JE.
    A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly.
    ACS Synthetic Biology 2015, 4(9), 975-986.  doi:10.1021/sb500366v
    Plasmids: Addgene kit #1000000061.

Two inputs, both from that kit:
    pYTK*.gb        one GenBank file per plasmid, gives the sequences
    YTK_Parts.xls   the parts list, gives the clean names and part types

Both live in reference/ytk_plasmids/, so this needs no arguments:
    python3 data/build_parts_table.py

Pass a directory to read them from somewhere else.

It writes data/ytk_parts.tsv next to this script.

Each pYTK plasmid carries one part between two BsaI sites. Cutting with BsaI
gives two pieces: the storage backbone and the part. The backbone is the piece
holding CamR, which no part carries -- two parts are longer than the backbone,
so picking by size would be wrong. Eight plasmids in the kit have no CamR at
all (see cut_out_part).

The table also holds four rows for the entry vector's own features
(ColE1, CamR, CamR Promoter, CamR Terminator) taken from pYTK001.gb. Those
four survive cloning and get labelled on finished plasmids.
"""
import sys
from pathlib import Path

from Bio import SeqIO
from Bio.Restriction import BsaI
from Bio.Seq import Seq
from pydna.dseqrecord import Dseqrecord
import pandas as pd

OUT = Path(__file__).resolve().parent / "ytk_parts.tsv"
PLASMIDS = Path(__file__).resolve().parents[1] / "reference" / "ytk_plasmids"

COLUMNS = ["name", "part_type", "feature_type", "strand",
           "junction_5", "junction_3", "sequence"]

# A circle has no natural first base, so a map has to pick one. SnapGene's own
# pYTK001 map starts here, and maps drawn by hand in SnapGene inherit that
# start. The entry vector is stored turned to this point, so plasmids built by
# this plugin line up with maps made by hand instead of looking different.
# These 40 bases sit inside ColE1, which stays in every finished plasmid.
ENTRY_VECTOR_START = "TCCTGTCGGGTTTCGCCACCTCTGACTTGAGCGTCGATTT"

# Entry vector features that stay in every assembled plasmid.
# Read from pYTK001.gb, which is correct. The project's own pYTK001.dna is not:
# it claims ColE1 spans the whole plasmid instead of 764 bp.
KEEP_FEATURES = {
    "ColE1": "rep_origin",
    "CamR": "CDS",
    "CamR Promoter": "promoter",
    "CamR Terminator": "terminator",
}


def feature_type_for(part_type):
    """Map a YTK part type to a SnapGene feature type."""
    t = str(part_type).strip().lower()
    if t == "2":
        return "promoter"
    if t.startswith("3"):
        return "CDS"
    if t.startswith("4"):
        return "terminator"
    return "misc_feature"


def revcomp(seq):
    return str(Seq(seq).reverse_complement()).upper()


def cut_out_part(plasmid_seq, cole1, camr):
    """Return (part sequence, 5' junction, 3' junction) for one pYTK plasmid.

    The part is the piece that does not hold CamR, the storage backbone's own
    marker. The type 8, 8a, 678 and cassette plasmids are kept on Amp, Kan or
    Spec instead and have no CamR anywhere: there the part itself is the
    bacterial origin plus marker, so it is the piece that holds ColE1.

    Neither landmark can be dropped. The GFP dropout is a part in pYTK047 and
    the leftover backbone in pYTK096, so looking at one piece alone can never
    tell the two apart -- what settles it is whether the plasmid has CamR.
    """
    frags = Dseqrecord(plasmid_seq, circular=True).cut(BsaI)
    if not frags:
        return None
    no_camr = [f for f in frags if camr not in str(f.seq).upper()]
    if len(no_camr) == 1:
        part = no_camr[0]
    else:
        part = [f for f in frags if cole1 in str(f.seq).upper()][0]
    j5 = part.seq.five_prime_end()[1].upper()
    j3 = revcomp(part.seq.three_prime_end()[1])
    return str(part.seq).upper(), j5, j3


def main(ytk_dir):
    ytk_dir = Path(ytk_dir).expanduser()
    parts_list = pd.read_excel(ytk_dir / "YTK_Parts.xls")

    entry = SeqIO.read(ytk_dir / "pYTK001.gb", "genbank")
    entry_seq = str(entry.seq).upper()

    def landmark(label):
        f = [f for f in entry.features
             if f.qualifiers.get("label", [""])[0] == label][0]
        return entry_seq[int(f.location.start):int(f.location.end)]

    cole1 = landmark("ColE1")
    camr = landmark("CamR")

    rows = []

    # The four entry vector features.
    for f in entry.features:
        label = f.qualifiers.get("label", [""])[0]
        if label not in KEEP_FEATURES:
            continue
        rows.append({
            "name": label,
            "part_type": "entry vector feature",
            "feature_type": KEEP_FEATURES[label],
            "strand": -1 if f.location.strand == -1 else 1,
            "junction_5": "",
            "junction_3": "",
            "sequence": entry_seq[int(f.location.start):int(f.location.end)],
        })

    # The whole entry vector, used as the cloning backbone. Turned so that it
    # begins at ENTRY_VECTOR_START. The feature slices above were taken before
    # this, while the GenBank coordinates still applied.
    start = entry_seq.index(ENTRY_VECTOR_START)
    rows.append({
        "name": "pYTK001", "part_type": "entry vector",
        "feature_type": "misc_feature", "strand": 1,
        "junction_5": "", "junction_3": "",
        "sequence": entry_seq[start:] + entry_seq[:start],
    })

    # One row per part.
    for _, r in parts_list.iterrows():
        plasmid = str(r["Plasmid Name"]).strip()
        gb = ytk_dir / f"{plasmid}.gb"
        if not gb.exists():
            print(f"skipped {plasmid}: no GenBank file", file=sys.stderr)
            continue
        cut = cut_out_part(str(SeqIO.read(gb, "genbank").seq).upper(),
                           cole1, camr)
        if cut is None:
            continue  # pYTK001 has no BsaI site; it is already in the table
        seq, j5, j3 = cut
        rows.append({
            "name": str(r["Part Description"]).strip(),
            "part_type": str(r["Part Type"]).strip(),
            "feature_type": feature_type_for(r["Part Type"]),
            "strand": 1,
            "junction_5": j5,
            "junction_3": j3,
            "sequence": seq,
        })

    pd.DataFrame(rows, columns=COLUMNS).to_csv(OUT, sep="\t", index=False)
    print(f"wrote {OUT} with {len(rows)} rows")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else PLASMIDS)
