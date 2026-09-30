"""The YTK flanks, shared by ytk-add-overhangs and ytk-design-primers.

Reading a flanked fragment outwards from the sequence:

    handle   BsmBI    BsaI      adapter   SEQUENCE   adapter    BsaI     BsmBI    handle
    actcgacaac CGTCTC atc GGTCTC a   T    ATG...TAA  ATCC    t GAGACC t GAGACG gttgtggtgt

BsaI cuts the part into an assembly and leaves the part-type overhangs, TATG
and ATCC for Type 3. BsmBI cuts the part out of the ordered fragment and leaves
TCGG and GACC, which is what the pYTK001 entry vector needs. The handles sit
outside both enzymes, so they are cut off and thrown away.

The two BsmBI sites have different spacers, `atc` on the left and `a` on the
right. That asymmetry is what makes the two overhangs come out different, so
the fragment can only go into the vector one way round. It is not a typo.

BsmBI outside and BsaI inside is how this plugin designs a fragment, and it is
the same for every part type. So ytk-clone's --enzyme picks which of those two
cuts is being simulated. It does not change what the flanks are.
"""
import csv
import sys
from pathlib import Path

DATA = Path(__file__).resolve().parents[1] / "data"
OVERHANG_TABLE = DATA / "ytk_overhangs.tsv"
PART_TYPE_TABLE = DATA / "ytk_part_types.tsv"

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

# The spare bases at each outer end, there only so BsmBI has room to cut near
# the end of a linear fragment. A PCR primer trims these; an ordered fragment
# keeps all ten.
LEFT_HANDLE = "actcgacaac"
RIGHT_HANDLE = "gttgtggtgt"

# Everything between a handle and the part-type adapter. FORWARD_SCAFFOLD is
# used by both the flanked fragment and the forward primer. RIGHT_SCAFFOLD is
# the reverse complement of REVERSE_SCAFFOLD, which a test checks.
FORWARD_SCAFFOLD = "CGTCTCatcGGTCTCa"
REVERSE_SCAFFOLD = "CGTCTCaGGTCTCa"
RIGHT_SCAFFOLD = "tGAGACCtGAGACG"

# The pad may be trimmed to make room for a longer binding region, but never
# below four bases or BsmBI loses its footing.
MIN_PAD = 4
FULL_PAD = len(LEFT_HANDLE)


def reverse_complement(sequence):
    return sequence.translate(COMPLEMENT)[::-1]


def _table(path):
    """Rows of a data table, keyed by part type.

    Comment rows are skipped. A row that stops early, because its last columns
    are empty, reads those columns as "" and not as None.
    """
    with open(path, newline="") as fh:
        rows = {}
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["part_type"].startswith("#"):
                continue
            rows[row["part_type"]] = {key: value or "" for key, value in row.items()}
        return rows


def adapters(part_type):
    """The adapter pair, and what the type holds, for one part type."""
    table = _table(OVERHANG_TABLE)
    if part_type not in table:
        sys.exit(f"unknown part type {part_type}. Known types: {', '.join(table)}")
    return table[part_type]


def junctions(part_type):
    """The overhang pair a part of this type must have.

    These come from the published part-type table, not from the adapters. The
    part type is defined by these four bases on each side and by nothing else,
    so reading them off a sequence would only ever agree with itself.
    """
    table = _table(PART_TYPE_TABLE)
    if part_type not in table:
        sys.exit(f"part type {part_type} has no published overhang pair. "
                 f"Known types: {', '.join(table)}")
    row = table[part_type]
    return row["upstream"].upper(), row["downstream"].upper()


def flank(sequence, adapters):
    """The full fragment to order, with both handles at their full length."""
    return (LEFT_HANDLE + FORWARD_SCAFFOLD + adapters["left_adapter"]
            + sequence
            + adapters["right_adapter"] + RIGHT_SCAFFOLD + RIGHT_HANDLE)


def insert_offset(adapters):
    """Where the sequence itself starts inside the flanked fragment."""
    return len(LEFT_HANDLE + FORWARD_SCAFFOLD + adapters["left_adapter"])


def forward_pad(length):
    """The forward pad, trimmed from its outer end."""
    return LEFT_HANDLE[-length:]


def reverse_pad(length):
    """The reverse pad, trimmed from its outer end.

    The spreadsheet only stores the top-strand handle, so the reverse pad is
    the reverse complement of it: acaccacaac.
    """
    return reverse_complement(RIGHT_HANDLE)[-length:]


def forward_primer(binding, adapters, pad=MIN_PAD):
    """A forward primer: pad, both enzyme sites, the left adapter, then the
    bases that stick to the template.

    For Type 3 the left adapter is a single T, and the sequence's own ATG
    completes the TATG overhang.
    """
    return forward_pad(pad) + FORWARD_SCAFFOLD + adapters["left_adapter"] + binding


def reverse_primer(binding, adapters, pad=MIN_PAD):
    """A reverse primer. Its adapter is the right adapter read on the other
    strand, so Type 3's ATCC appears as GGAT.
    """
    return (reverse_pad(pad) + REVERSE_SCAFFOLD
            + reverse_complement(adapters["right_adapter"]) + binding)


def constant_length(adapters, pad=MIN_PAD):
    """How many bases of a primer are scaffold rather than binding region."""
    return (len(forward_primer("", adapters, pad)),
            len(reverse_primer("", adapters, pad)))
