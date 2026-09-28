"""Read and write plasmid maps, and compare circular sequences.

Shared by more than one skill, so it lives at the repo root.

Two formats. SnapGene .dna is written here by hand, byte by byte, and needs
nothing installed. GenBank goes through Biopython. write_map and read_map pick
between them from the file suffix, so a caller names a file and nothing else.

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

The .dna half imports nothing, so it always works. Only the GenBank functions
need Biopython, and they import it when called.
"""
import csv
import datetime
import struct
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import deps

CIRCULAR = 0x01
DOUBLE_STRANDED = 0x02


def require(package, module=None):
    """Stop with the pip command when a package is missing.

    Nothing is installed here, because pip would fetch the newest release
    rather than the version pinned in requirements.txt.
    """
    deps.require(module or package)


# --- reading the gene file ---------------------------------------------

NAME_HEADERS = {"name", "gene", "id", "gene_name", "gene name"}
SEQUENCE_HEADERS = {"sequence", "seq", "dna"}
PLASMID_HEADERS = {"plasmid", "plasmid_name", "plasmid name", "construct"}
NOTES_HEADERS = {"notes", "note"}
PART_TYPE_HEADERS = {"part_type", "part type", "type"}


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


def _rows(path):
    """Every non-blank row, and the header if there is one.

    utf-8-sig, because a CSV saved by Excel starts with a byte order mark.
    Without it the first header reads as "﻿name", matches nothing, and the
    header row is then parsed as a gene. Both readers go through here, so they
    agree on what a header is and neither can drift from the other.
    """
    with open(path, newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.reader(fh) if any(cell.strip() for cell in r)]

    if not rows:
        return [], None

    header = [cell.strip().lower() for cell in rows[0]]
    if set(header) & NAME_HEADERS and set(header) & SEQUENCE_HEADERS:
        return rows[1:], header
    return rows, None


def _column_at(header, wanted):
    return next((i for i, h in enumerate(header) if h in wanted), None)


def _read_gene_csv(path):
    rows, header = _rows(path)
    if not rows:
        return []

    if header:
        name_at = _column_at(header, NAME_HEADERS)
        sequence_at = _column_at(header, SEQUENCE_HEADERS)
        plasmid_at = _column_at(header, PLASMID_HEADERS)
    else:
        name_at, sequence_at, plasmid_at = 0, 1, 2

    genes = []
    for row in rows:
        plasmid = None
        if plasmid_at is not None and len(row) > plasmid_at:
            plasmid = row[plasmid_at].strip() or None
        genes.append((row[name_at].strip(), row[sequence_at].strip().upper(), plasmid))
    return genes


def read_column(path, wanted):
    """An optional extra column of a gene CSV, as {gene name: value}.

    Used for the `notes` and `part_type` columns that ytk-add-overhangs writes
    and later steps carry through. Returns {} for FASTA, GenBank, .dna, a CSV
    with no header, or a CSV without that column - the columns are optional
    everywhere, so a missing one is not an error.

    The value rides on the same row as the sequence it describes, which is the
    point: there is no join key, so a note can never end up attached to a
    different sequence than the one it was written for.
    """
    path = Path(path)
    if path.suffix.lower() != ".csv":
        return {}

    rows, header = _rows(path)
    if not header:
        return {}
    name_at = _column_at(header, NAME_HEADERS)
    value_at = _column_at(header, wanted)
    if value_at is None:
        return {}

    found = {}
    for row in rows:
        if len(row) > value_at:
            found[row[name_at].strip()] = row[value_at].strip()
    return found


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


# --- GenBank -----------------------------------------------------------

GENBANK_SUFFIXES = {".gb", ".gbk"}


def _genbank_feature(feature, length):
    from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature

    strand = feature.get("strand", 1)
    wrap_end = feature.get("wrap_end")
    if wrap_end:
        # A feature that runs past the end of the circle is a join() of the
        # two pieces. Biopython wants the parts in 5' to 3' order, which on
        # the reverse strand is the other way round.
        parts = [FeatureLocation(feature["start"], length, strand),
                 FeatureLocation(0, wrap_end, strand)]
        if strand == -1:
            parts.reverse()
        location = CompoundLocation(parts)
    else:
        location = FeatureLocation(feature["start"], feature["end"], strand)
    return SeqFeature(location, type=feature["type"],
                      qualifiers={"label": [feature["name"]]})


def write_genbank(path, sequence, circular, notes_type="Natural", description="",
                  features=None, created_by="IOCB Prague"):
    """Write a GenBank file, taking the same feature dicts as write_dna.

    notes_type is SnapGene's own field and has no GenBank equivalent, so it is
    accepted and ignored. That keeps the signature the same as write_dna, which
    is what lets write_map forward to either one.
    """
    require("biopython", "Bio")
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    name = Path(path).stem
    record = SeqRecord(Seq(sequence.upper()), id=name, name=name,
                       description=description)
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "circular" if circular else "linear"
    record.annotations["date"] = datetime.date.today().strftime("%d-%b-%Y").upper()
    record.annotations["source"] = created_by
    record.features = [_genbank_feature(f, len(sequence)) for f in features or []]
    SeqIO.write(record, str(path), "genbank")


def read_genbank(path):
    """Read a GenBank file into the same dict shape read_dna returns.

    double_stranded is always True: a GenBank record has no such flag, because
    the format assumes double-stranded DNA. Anything reporting that check has
    to say it did not really run.
    """
    require("biopython", "Bio")
    from Bio import SeqIO

    record = SeqIO.read(str(path), "genbank")
    return {"sequence": str(record.seq).upper(),
            "circular": record.annotations.get("topology") == "circular",
            "double_stranded": True}


# --- either format -----------------------------------------------------


def _is_genbank(path):
    return Path(path).suffix.lower() in GENBANK_SUFFIXES


def write_map(path, sequence, circular, **kwargs):
    """Write a map as GenBank or SnapGene, picked from the file suffix.

    Reads the file back and compares the sequence. A writer that dropped or
    reordered bases would otherwise leave a file that opens perfectly and is
    wrong, and nothing downstream would catch it.
    """
    writer = write_genbank if _is_genbank(path) else write_dna
    writer(path, sequence, circular, **kwargs)

    written = read_map(path)["sequence"]
    if written != sequence.upper():
        raise ValueError(
            f"{path}: the file read back different from what was written "
            f"({len(written)} bp against {len(sequence)}). The file is wrong. "
            f"Send this to whoever maintains the plugin.")


def read_map(path):
    """Read a map from GenBank or SnapGene, picked from the file suffix."""
    return read_genbank(path) if _is_genbank(path) else read_dna(path)


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
