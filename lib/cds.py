"""Does this sequence fit the part type asked for?

Every check here is driven by the `coding` column of data/ytk_overhangs.tsv,
which says how much of a protein the type holds:

    full    a whole coding sequence: ATG at the start, a stop codon at the end
    start   the first half of a protein fusion: ATG, but no stop codon
    end     the second half: mid-protein at both ends
    ""      not a coding sequence at all

For a type with no coding rule there is nothing to check beyond the two things
that are true of any DNA: it is not empty, and its letters are A, C, G and T.
A promoter has no testable shape and no reading frame, so nothing is invented
for it.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import enzymes

STOP_CODONS = ("TAA", "TAG", "TGA")

# One fixed wording, so add_overhangs can pick this fault out of the list when
# --no-stop-codon says the missing stop is on purpose.
NO_STOP_CODON = "does not end with a stop codon"


def problems(sequence, coding):
    """Every reason this does not fit the type. Empty means it does.

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
    if coding and len(sequence) % 3:
        found.append(f"is {len(sequence)} bases, not a whole number of codons")
    # A left adapter of a single T relies on the gene's own ATG to complete the
    # TATG overhang, so without that ATG the junction is simply wrong.
    if coding in ("full", "start") and not sequence.startswith("ATG"):
        found.append(f"starts with {sequence[:3]}, not ATG")
    ends_in_stop = sequence[-3:] in STOP_CODONS
    if coding in ("full", "end") and not ends_in_stop:
        found.append(NO_STOP_CODON)
    if coding == "start" and ends_in_stop:
        found.append("is the first half of a protein fusion, so it must not end "
                     "with a stop codon")
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


def warnings(sequence, coding):
    """Things worth saying about a sequence that does fit its type."""
    found = [f"holds a {enzyme} site, which breaks Golden Gate"
             for enzyme in enzymes.in_sequence(sequence)]
    if not coding:
        return found
    position = internal_stop(sequence)
    if position is not None:
        found.append(f"has a stop codon at base {position + 1}, before the end")
    return found
