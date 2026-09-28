"""Read sequences from any of the four formats the lab uses.

FASTA and CSV hold one sequence per record, so there is nothing to choose.
GenBank and SnapGene files hold a whole map with many features, so the coding
sequence has to be picked out. That choice is never guessed: when it is not
obvious the run stops and asks.
"""
import csv
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import deps
import snapgene as sg

DATA = Path(__file__).resolve().parents[1] / "data"
PARTS_TABLE = DATA / "ytk_parts.tsv"
BACKBONE_TABLE = DATA / "backbone_features.tsv"

PLAIN_FORMATS = {".fa", ".fasta", ".csv"}
GENBANK_FORMATS = {".gb", ".gbk"}

# Feature types that are never the coding sequence, so they are not worth
# offering.
NOT_A_GENE = {"source", "primer_bind", "rep_origin", "terminator", "promoter",
              "RBS", "protein_bind"}

# Below this, a sequence match against the parts table is not worth trusting.
SHORTEST_PART_MATCH = 60


class AmbiguousCDS(Exception):
    """The file holds several possible coding sequences, or none.

    Carries the candidates so a skill can list them and ask.
    """

    def __init__(self, path, candidates):
        self.path = path
        self.candidates = candidates
        super().__init__(f"{path.name}: cannot tell which feature is the "
                         f"coding sequence")

    def table(self):
        lines = [f"{'name':28s} {'type':14s} {'length':>7s} {'position':>14s}"]
        for feature in self.candidates:
            lines.append(f"{feature['name'][:28]:28s} {feature['type'][:14]:14s} "
                         f"{feature['end'] - feature['start']:7d} "
                         f"{feature['start'] + 1:6d}-{feature['end']:<7d}")
        return "\n".join(lines)


def read(path, feature=None):
    """Sequences as (name, sequence, plasmid name or None)."""
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix in PLAIN_FORMATS:
        return sg.read_genes(path)
    if suffix in GENBANK_FORMATS:
        return _read_genbank(path, feature)
    if suffix == ".dna":
        return [_read_snapgene(path, feature)]
    sys.exit(f"{path.name}: expected .fa, .fasta, .csv, .gb, .gbk or .dna")


def _read_genbank(path, feature):
    deps.require("Bio")
    from Bio import SeqIO

    found = []
    for record in SeqIO.parse(str(path), "genbank"):
        features = [
            {"name": _genbank_label(item) or record.id,
             "type": item.type,
             "start": int(item.location.start),
             "end": int(item.location.end),
             "strand": item.location.strand or 1}
            for item in record.features
        ]
        whole = str(record.seq).upper()
        chosen = _choose(path, _mark_known_parts(features, whole), feature)
        sequence = whole[chosen["start"]:chosen["end"]]
        if chosen["strand"] == -1:
            sequence = _reverse_complement(sequence)
        found.append((chosen["name"], sequence, None))
    return found


def _genbank_label(item):
    for key in ("label", "gene", "product"):
        if key in item.qualifiers:
            return item.qualifiers[key][0]
    return None


def _read_snapgene(path, feature):
    contents = sg.read_dna(path)
    whole = contents["sequence"].upper()
    features = _snapgene_features(contents.get("features_xml", ""))
    chosen = _choose(path, _mark_known_parts(features, whole), feature)
    sequence = whole[chosen["start"]:chosen["end"]]
    if chosen["strand"] == -1:
        sequence = _reverse_complement(sequence)
    return (chosen["name"], sequence, None)


def _snapgene_features(xml):
    """Features from the XML that SnapGene keeps in chunk 10.

    A feature carries its label in a <Q name="label"> block when it has one,
    and its span in one or more <Segment range="start-end"> children, counted
    from 1 and including both ends.
    """
    if not xml:
        return []
    found = []
    for item in ET.fromstring(xml).findall("Feature"):
        spans = [segment.get("range", "") for segment in item.findall("Segment")]
        edges = [int(part) for span in spans for part in span.split("-") if part]
        if not edges:
            continue
        label = item.find("Q[@name='label']/V")
        found.append({
            "name": (label.get("text") if label is not None else None)
                    or item.get("name", "unnamed"),
            "type": item.get("type", "misc_feature"),
            "start": min(edges) - 1,
            "end": max(edges),
            "strand": -1 if item.get("directionality") == "2" else 1,
        })
    return found


def _choose(path, features, wanted):
    """Which feature holds the coding sequence.

    A name given on the command line wins. Otherwise a single coding feature is
    taken. Anything less certain stops the run, because picking the wrong
    feature would put the wrong DNA in a plasmid and nothing downstream would
    notice.
    """
    if wanted:
        matches = [f for f in features if f["name"].lower() == wanted.lower()]
        if not matches:
            sys.exit(f"{path.name}: no feature called {wanted}. Found: "
                     + ", ".join(sorted(f["name"] for f in features)))
        return matches[0]

    # Preferring the one feature typed CDS looks reasonable and is not safe: in
    # a real part plasmid the only CDS is the CamR marker, while the insert is
    # labelled misc_feature. So known backbone machinery is set aside first.
    possible = sorted((f for f in features if f["type"] not in NOT_A_GENE),
                      key=lambda f: f["start"] - f["end"])
    candidates = [f for f in possible if not f.get("known_part")]
    if len(candidates) == 1:
        return candidates[0]
    # Everything recognised means there is no insert here, so offer the lot.
    raise AmbiguousCDS(path, candidates or possible)


def _load_parts():
    """Sequences of the published YTK parts, both strands, longest first."""
    known = {}
    with open(PARTS_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sequence = row["sequence"].upper()
            if len(sequence) >= SHORTEST_PART_MATCH:
                known[sequence] = row["name"]
                known[_reverse_complement(sequence)] = row["name"]
    return known


def _load_backbone_names():
    """(normalised name, exact?) for features that belong to a vector."""
    names = []
    with open(BACKBONE_TABLE, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["name"].startswith("#") or not row["name"]:
                continue
            names.append((_normalise(row["name"]), row["match"] == "exact"))
    return names


def _normalise(name):
    return "".join(character for character in name.lower()
                   if character.isalnum())


def _mark_known_parts(features, sequence):
    """Flag the features that are vector machinery rather than an insert.

    A sequence match against the published parts table is the strong test,
    because it holds whatever the feature happens to be labelled. The name
    catalogue is the backstop for vectors that are not from the toolkit.
    """
    parts = _load_parts()
    catalogue = _load_backbone_names()
    for feature in features:
        span = sequence[feature["start"]:feature["end"]].upper()
        name = _normalise(feature["name"])
        by_sequence = parts.get(span)
        by_name = next((entry for entry, exact in catalogue
                        if name == entry or (not exact and entry in name)), None)
        if by_sequence or by_name:
            feature["known_part"] = by_sequence or by_name
    return features


def _reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]
