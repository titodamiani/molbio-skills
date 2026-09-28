"""Put labels on a SnapGene plasmid map.

Reads a .dna file, looks for known YTK parts in it, and writes the file back
with those parts labelled.

    python3 annotate.py out/pTP412.dna
    python3 annotate.py out/*.dna --genes genes.fasta

Parts come from the shared table data/ytk_parts.tsv, which was built from the
published GenBank files. Labels are never copied from another map file,
because those can be wrong.

With --genes, the cloned gene is labelled too. Each gene in that file is
searched for, and the one that is present gets its own name on the map.

Parts are found by searching for their sequence, so the labels stay right no
matter where the map starts.
"""
import argparse
import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import sequences  # noqa: E402
import snapgene as sg

PARTS_TABLE = ROOT / "data" / "ytk_parts.tsv"

COLORS = {
    "CDS": "#993366",
    "promoter": "#ffffff",
    "terminator": "#a6acb3",
    "rep_origin": "#ffff00",
    "misc_feature": "#c6c9d1",
}
GENE_COLOR = "#66ccff"


def load_parts():
    with open(PARTS_TABLE, newline="") as fh:
        return [r for r in csv.DictReader(fh, delimiter="\t")
                if r["part_type"] != "entry vector"]


def make_feature(name, feature_type, start, length, strand, color, total):
    """Build one feature, marking it if it runs past the end of the circle."""
    end = start + length
    return {
        "name": name,
        "type": feature_type,
        "start": start % total,
        "end": end if end <= total else total,
        "strand": strand,
        "color": color,
        "wrap_end": (end - total) if end > total else 0,
    }


def build_features(plasmid, parts, genes):
    total = len(plasmid)
    features = []

    for name, gene in genes:
        at = sg.find_in_circle(plasmid, gene)
        if at is not None:
            features.append(make_feature(name, "CDS", at, len(gene), 1,
                                         GENE_COLOR, total))

    for part in parts:
        at = sg.find_in_circle(plasmid, part["sequence"])
        if at is None:
            continue
        feature_type = part["feature_type"]
        features.append(make_feature(
            part["name"], feature_type, at, len(part["sequence"]),
            int(part["strand"]), COLORS.get(feature_type, "#c6c9d1"), total))

    # Cut sites and fusion scars are left off on purpose. SnapGene shows those
    # live under Enzymes and Common Features, so a fixed label would go stale.
    features.sort(key=lambda f: f["start"])
    return features


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("plasmids", nargs="+", help=".dna files to label")
    ap.add_argument("--genes", help="FASTA or CSV file of the cloned genes")
    args = ap.parse_args()

    parts = load_parts()
    genes = [(n, s) for n, s, _ in sequences.read(args.genes)] if args.genes else []

    for path in args.plasmids:
        record = sg.read_dna(path)
        features = build_features(record["sequence"], parts, genes)
        sg.write_dna(path, record["sequence"], circular=record["circular"],
                     notes_type="Synthetic" if record["circular"] else "Natural",
                     description="synthetic circular DNA." if record["circular"] else "",
                     features=features)
        names = ", ".join(f["name"] for f in features)
        print(f"{Path(path).name}: {len(features)} labels -- {names}")


if __name__ == "__main__":
    main()
