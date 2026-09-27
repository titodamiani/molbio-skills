"""Is this a YTK Type 3 coding sequence?

Type 3 holds a whole coding sequence: it starts at ATG, ends with a stop codon,
and is a whole number of codons. The plugin only builds Type 3, so there is
nothing here about the other part types.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import enzymes

STOP_CODONS = ("TAA", "TAG", "TGA")


def problems(sequence):
    """Every reason this is not a Type 3 coding sequence. Empty means it is.

    All of them, not just the first: someone fixing a sequence wants the whole
    list in one go.
    """
    sequence = sequence.upper()
    found = []
    if not sequence:
        return ["is empty"]
    odd = sorted(set(sequence) - set("ACGT"))
    if odd:
        found.append(f"holds {', '.join(odd)} as well as A, C, G and T")
    if not sequence.startswith("ATG"):
        found.append(f"starts with {sequence[:3]}, not ATG")
    if sequence[-3:] not in STOP_CODONS:
        found.append(f"ends with {sequence[-3:]}, which is not a stop codon")
    if len(sequence) % 3:
        found.append(f"is {len(sequence)} bases, not a whole number of codons")
    return found


def internal_stop(sequence):
    """Position of the first stop codon before the end, or None.

    A stop in the middle means the reading frame is wrong, or the sequence is
    not the coding sequence someone thought it was.
    """
    sequence = sequence.upper()
    for start in range(0, len(sequence) - 3, 3):
        if sequence[start:start + 3] in STOP_CODONS:
            return start
    return None


def warnings(sequence):
    """Things worth saying about a sequence that is otherwise a valid Type 3."""
    found = [f"holds a {enzyme} site, which breaks Golden Gate"
             for enzyme in enzymes.in_sequence(sequence)]
    position = internal_stop(sequence)
    if position is not None:
        found.append(f"has a stop codon at base {position + 1}, before the end")
    return found
