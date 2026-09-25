"""Build data/ytk_parts.tsv from the published MoClo Yeast Toolkit files.

Source (supplementary material of):
    Lee ME, DeLoache WC, Cervantes B, Dueber JE.
    A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly.
    ACS Synthetic Biology 2015, 4(9), 975-986.  doi:10.1021/sb500366v
    Plasmids: Addgene kit #1000000061.

Two inputs, both from that kit:
    pYTK*.gb        one GenBank file per plasmid, gives the sequences
    YTK_Parts.xls   the parts list, gives the clean names and part types

Run it like this:
    python3 data/build_parts_table.py ~/Downloads/ytk/plasmids

It writes data/ytk_parts.tsv next to this script.

Each pYTK plasmid carries one part between two BsaI sites. Cutting with BsaI
gives two pieces: the entry vector (it holds ColE1) and the part. The part is
the piece without ColE1 -- two parts are longer than the vector, so picking by
size would be wrong.

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


def cut_out_part(plasmid_seq, cole1):
    """Return (part sequence, 5' junction, 3' junction) for one pYTK plasmid."""
    frags = Dseqrecord(plasmid_seq, circular=True).cut(BsaI)
    if not frags:
        return None
    part = [f for f in frags if cole1 not in str(f.seq).upper()][0]
    j5 = part.seq.five_prime_end()[1].upper()
    j3 = revcomp(part.seq.three_prime_end()[1])
    return str(part.seq).upper(), j5, j3


def main(ytk_dir):
    ytk_dir = Path(ytk_dir).expanduser()
    parts_list = pd.read_excel(ytk_dir / "YTK_Parts.xls")

    entry = SeqIO.read(ytk_dir / "pYTK001.gb", "genbank")
    entry_seq = str(entry.seq).upper()
    cole1 = [entry_seq[int(f.location.start):int(f.location.end)]
             for f in entry.features
             if f.qualifiers.get("label", [""])[0] == "ColE1"][0]

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
        cut = cut_out_part(str(SeqIO.read(gb, "genbank").seq).upper(), cole1)
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
    main(sys.argv[1] if len(sys.argv) > 1 else "~/Downloads/ytk/plasmids")
