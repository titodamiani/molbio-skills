"""Read and write SnapGene .dna files, and compare circular sequences.

Shared by more than one skill, so it lives at the repo root.

A .dna file is a chain of chunks. Each chunk is 1 byte for the chunk type,
then 4 bytes (big-endian) for the payload length, then the payload:

    9   header      the "SnapGene" cookie and version numbers
    0   sequence    1 flags byte, then the bases as ASCII
    6   notes       XML
    10  features    XML

Flags byte on chunk 0: bit 0x01 means circular, bit 0x02 means double
stranded. With 0x02 off, SnapGene treats the file as a single strand and
hides the Enzymes tab. So linear DNA is 0x02 and a plasmid is 0x03.

Feature positions in the XML are 1-based and include both ends. A feature
that runs past the end of a circle is written with its end number lower than
its start. directionality="2" means the reverse strand.

Nothing here imports a third-party package, so it always works.
"""
import csv
import datetime
import struct
import subprocess
import sys
from pathlib import Path

CIRCULAR = 0x01
DOUBLE_STRANDED = 0x02


def require(package, module=None):
    """Import a package, installing it with pip the first time if it is absent."""
    module = module or package
    try:
        __import__(module)
    except ImportError:
        print(f"installing {package} ...", file=sys.stderr)
        subprocess.check_call([sys.executable, "-m", "pip", "install", package])
        __import__(module)


# --- reading the gene file ---------------------------------------------

NAME_HEADERS = {"name", "gene", "id", "gene_name", "gene name"}
SEQUENCE_HEADERS = {"sequence", "seq", "dna"}
PLASMID_HEADERS = {"plasmid", "plasmid_name", "plasmid name", "construct"}


def read_genes(path):
    """Read the gene file. Returns a list of (name, sequence, plasmid name).

    The plasmid name is None when the file does not give one.

    Two kinds of file are accepted, and nothing else:

        FASTA   .fa or .fasta, which cannot carry a plasmid name
        CSV     .csv, with a name column, a sequence column, and an optional
                plasmid name column

    A CSV may have a header row or not, and capitals in the header do not
    matter. Every skill reads gene files through here, so they all agree on
    what a gene file is.
    """
    path = Path(path)
    suffix = path.suffix.lower()

    if suffix in (".fa", ".fasta"):
        require("biopython", "Bio")
        from Bio import SeqIO
        genes = [(r.id, str(r.seq).upper(), None) for r in SeqIO.parse(path, "fasta")]
    elif suffix == ".csv":
        genes = _read_gene_csv(path)
    else:
        sys.exit(f"{path.name}: only .fa, .fasta and .csv files are accepted, "
                 f"not {suffix or 'a file without an extension'}")

    _check_genes(path, genes)
    return genes


def _read_gene_csv(path):
    with open(path, newline="") as fh:
        rows = [r for r in csv.reader(fh) if any(cell.strip() for cell in r)]

    header = [cell.strip().lower() for cell in rows[0]]
    name_at, sequence_at, plasmid_at = 0, 1, 2
    if set(header) & NAME_HEADERS and set(header) & SEQUENCE_HEADERS:
        name_at = next(i for i, h in enumerate(header) if h in NAME_HEADERS)
        sequence_at = next(i for i, h in enumerate(header) if h in SEQUENCE_HEADERS)
        plasmid_at = next((i for i, h in enumerate(header) if h in PLASMID_HEADERS), None)
        rows = rows[1:]

    genes = []
    for row in rows:
        plasmid = None
        if plasmid_at is not None and len(row) > plasmid_at:
            plasmid = row[plasmid_at].strip() or None
        genes.append((row[name_at].strip(), row[sequence_at].strip().upper(), plasmid))
    return genes


def _check_genes(path, genes):
    """Stop on anything that would give a wrong or overwritten file."""
    if not genes:
        sys.exit(f"{path.name}: no sequences found")

    for name, sequence, _ in genes:
        odd = sorted(set(sequence) - set("ACGT"))
        if odd:
            sys.exit(f"{name}: the sequence contains {', '.join(odd)}. "
                     f"Only A, C, G and T are allowed. Nothing was written.")

    # Two rows sharing a name would write one file over the other, and the
    # loss would be silent. Better to stop and let someone fix the list.
    for label, values in (("gene name", [g[0] for g in genes]),
                          ("plasmid name", [g[2] for g in genes if g[2]])):
        repeated = sorted({v for v in values if values.count(v) > 1})
        if repeated:
            sys.exit(f"{path.name}: the {label} {', '.join(repeated)} is used "
                     f"more than once. Each one must be different, or one file "
                     f"would overwrite another. Nothing was written.")


# --- writing -----------------------------------------------------------


def _chunk(chunk_type, payload):
    return struct.pack(">B", chunk_type) + struct.pack(">I", len(payload)) + payload


def _header_chunk():
    # sequence type 1 = DNA, export version 15, import version 20
    return _chunk(9, struct.pack(">8sHHH", b"SnapGene", 1, 15, 20))


def _sequence_chunk(sequence, circular):
    flags = DOUBLE_STRANDED | (CIRCULAR if circular else 0)
    return _chunk(0, struct.pack(">B", flags) + sequence.upper().encode("ascii"))


def _notes_chunk(seq_type, description, created_by):
    today = datetime.date.today().strftime("%Y.%m.%d")
    xml = (
        "<Notes>\n"
        f"<Type>{seq_type}</Type>\n"
        f"<Description>{description}</Description>\n"
        f'<Created UTC="0:0:0">{today}</Created>\n'
        f'<LastModified UTC="0:0:0">{today}</LastModified>\n'
        f"<CreatedBy>{created_by}</CreatedBy>\n"
        "<SequenceClass>UNA</SequenceClass>\n"
        "<TransformedInto>unspecified</TransformedInto>\n"
        "</Notes>\n"
    )
    return _chunk(6, xml.encode("utf-8"))


def _feature_xml(feature_id, feature):
    reverse = feature.get("strand", 1) == -1
    directionality = ' directionality="2"' if reverse else ""
    swapped = ' swappedSegmentNumbering="1"' if reverse else ""
    start = feature["start"] + 1
    # A feature that runs past the end of a circle gets an end number lower
    # than its start. SnapGene reads that as a wrap.
    end = feature["wrap_end"] if feature.get("wrap_end") else feature["end"]
    return (
        f'<Feature recentID="{feature_id}" name="{feature["name"]}"{directionality} '
        f'type="{feature["type"]}" allowSegmentOverlaps="0" '
        f'consecutiveTranslationNumbering="1"{swapped}>'
        f'<Segment range="{start}-{end}" color="{feature.get("color", "#a6acb3")}" '
        f'type="standard"/>'
        f'<Q name="label"><V text="{feature["name"]}"/></Q>'
        "</Feature>"
    )


def _features_chunk(features):
    body = "".join(_feature_xml(i, f) for i, f in enumerate(features))
    xml = f'<?xml version="1.0"?><Features nextValidID="{len(features)}">{body}</Features>'
    return _chunk(10, xml.encode("utf-8"))


def write_dna(path, sequence, circular, notes_type="Natural", description="",
              features=None, created_by="IOCB Prague"):
    """Write a SnapGene .dna file.

    features is a list of dicts with keys name, type, start (0-based),
    end (exclusive), strand, color, and wrap_end for a feature that runs past
    the end of a circle.
    """
    with open(path, "wb") as fh:
        fh.write(_header_chunk())
        fh.write(_sequence_chunk(sequence, circular))
        fh.write(_notes_chunk(notes_type, description, created_by))
        if features:
            fh.write(_features_chunk(features))


# --- reading -----------------------------------------------------------


def read_dna(path):
    """Read a .dna file. Returns a dict with sequence, circular, double_stranded."""
    with open(path, "rb") as fh:
        data = fh.read()

    result = {}
    at = 0
    while at + 5 <= len(data):
        chunk_type = data[at]
        length = struct.unpack(">I", data[at + 1:at + 5])[0]
        payload = data[at + 5:at + 5 + length]
        at += 5 + length
        if chunk_type == 0:
            flags = payload[0]
            result["sequence"] = payload[1:].decode("ascii").upper()
            result["circular"] = bool(flags & CIRCULAR)
            result["double_stranded"] = bool(flags & DOUBLE_STRANDED)
        elif chunk_type == 10:
            result["features_xml"] = payload.decode("utf-8")
    return result


# --- circles -----------------------------------------------------------


def find_in_circle(circle, fragment):
    """Return where fragment starts in a circular sequence, or None.

    Raises if the fragment is there more than once, because then any answer
    would be a guess.
    """
    circle = circle.upper()
    fragment = fragment.upper()
    doubled = circle + circle
    hits = []
    at = doubled.find(fragment)
    while at != -1 and at < len(circle):
        hits.append(at)
        at = doubled.find(fragment, at + 1)
    if not hits:
        return None
    if len(hits) > 1:
        raise ValueError(f"found {len(hits)} times, expected once")
    return hits[0]


def rotate_to(sequence, landmark):
    """Turn a circular sequence so that the landmark becomes base 1."""
    offset = find_in_circle(sequence, landmark)
    if offset is None:
        raise ValueError("landmark not found, cannot turn the circle")
    return sequence[offset:] + sequence[:offset]


def rotation_of(reference, other):
    """How far other is turned relative to reference, or None if they differ.

    Two files can describe the same plasmid while starting at different points
    on the circle. Comparing them as plain text would wrongly report a
    difference, so compare them as circles instead.
    """
    if len(reference) != len(other):
        return None
    doubled = other.upper() + other.upper()
    at = doubled.find(reference.upper())
    return at if 0 <= at < len(other) else None


def same_circle(a, b):
    """True when two sequences describe the same circular DNA."""
    return rotation_of(a, b) is not None
