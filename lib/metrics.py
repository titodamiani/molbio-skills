"""Measure a sequence. Nothing here changes one, ever.

Everything ytk-cds-qc can say about a sequence without touching it, in one
place, so the numbers are written into input.csv once and every later step reads
them off the row instead of measuring again or reading a sentence apart.

The numbers are here to decide which fixing skill a gene needs:

    a BsmBI or BsaI count above zero    ->  ytk-remove-cut-sites
    a low cai_scer                      ->  ytk-codon-optimise, perhaps
    a long run or a long repeat         ->  nothing here fixes it; say so

There is one CAI column and it is against yeast. CAI cannot say which organism a
gene was optimised for, however much it looks like it should: the index is a
geometric mean of each codon's share of its own family, so a table with more even
codon usage gives every sequence a higher score. The human table is flatter than
the yeast one (mean weight 0.73 against 0.66), and a random ORF that nothing
optimised scores 0.64 against human and 0.60 against yeast - the same ranking the
real genes in tests/data/ get. Columns for other organisms would have looked like
an answer and been an artefact of the tables.

GC and repeats are measured and reported and never fixed. Fixing either needs a
search over the whole sequence and a target nobody has set, and a number on a
table lets a person decide for themselves.

CAI comes from Biopython's Bio.SeqUtils.CodonAdaptationIndex, which implements
Sharp & Li (Nucleic Acids Research 15(3): 1281-1295 (1987)) and skips ATG, TGG
and the stop codon, as the published index does. Writing that by hand is a good
way to end up with a number nobody else can reproduce.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import deps
import enzymes

deps.require("Bio", "python_codon_tables")
import python_codon_tables as pct           # noqa: E402
from Bio.SeqUtils import CodonAdaptationIndex  # noqa: E402

# The shortest exact repeat worth reporting. Synthesis companies start asking
# questions around here.
REPEAT = 15
GC_WINDOW = 50

YEAST = 4932

_INDEX = {}


def gc(sequence):
    """Percent G plus C, or None for an empty sequence."""
    if not sequence:
        return None
    return round(100.0 * sum(sequence.count(base) for base in "GC")
                 / len(sequence), 1)


def gc_windows(sequence):
    """The lowest and highest GC of any GC_WINDOW-base stretch.

    A gene can sit at a comfortable 45% overall and still hold a 50-base window
    at 20%, which is the part a synthesis company struggles with. The whole-
    sequence figure hides exactly that, so both are reported.
    """
    if len(sequence) < GC_WINDOW:
        return gc(sequence), gc(sequence)
    windows = [gc(sequence[at:at + GC_WINDOW])
               for at in range(len(sequence) - GC_WINDOW + 1)]
    return min(windows), max(windows)


def longest_run(sequence):
    """The longest run of one base, as (base, length, 1-based position)."""
    best = ("", 0, 0)
    start = 0
    for at in range(1, len(sequence) + 1):
        if at == len(sequence) or sequence[at] != sequence[start]:
            if at - start > best[1]:
                best = (sequence[start], at - start, start + 1)
            start = at
    return best


def longest_repeat(sequence):
    """The longest exact repeat, as (length, first position, second position).

    (0, None, None) when nothing repeats over REPEAT bases. Every REPEAT-base
    window goes in a dict; the first window to come round twice is extended
    while the bases still agree. Any repeat longer than REPEAT contains such a
    window, so this finds the longest one.

    Measured on every sequence because the argmax rewrite in lib/codons.py
    creates repeats rather than removing them: the same peptide motif twice
    becomes the same DNA twice. Nothing here fixes that. Saying so is the point.
    """
    seen = {}
    best = (0, None, None)
    for at in range(len(sequence) - REPEAT + 1):
        window = sequence[at:at + REPEAT]
        if window not in seen:
            seen[window] = at
            continue
        first = seen[window]
        length = REPEAT
        while (at + length < len(sequence)
               and sequence[first + length] == sequence[at + length]):
            length += 1
        if length > best[0]:
            best = (length, first + 1, at + 1)
    return best


def _index(taxid):
    """Biopython's CAI index, built from a usage table rather than from genes.

    CodonAdaptationIndex counts the codons in a set of reference genes and works
    out each codon's share of its own family. python_codon_tables gives those
    shares directly, so the index is built on one throwaway codon - which sets
    up the genetic code it needs - and then every weight is written in. It is a
    dict subclass, so clear() and update() are its own public API and not a poke
    at insides.

    The three stop codons are left out, and that is what makes the number the
    published one. calculate() skips a stop codon only when the index has no
    weight for it; an index built the usual way, from whole genes, does carry
    stop weights, so it scores the terminator as though a gene could have chosen
    a different one. Sharp & Li exclude it. With the stops left out here, a
    sequence of nothing but top codons scores exactly 1.0, which is the check
    that this is right.
    """
    if taxid not in _INDEX:
        index = CodonAdaptationIndex(["ATG"])
        index.clear()
        index.update({codon: share / max(family.values())
                      for amino, family in pct.get_codons_table(taxid).items()
                      if amino != "*"
                      for codon, share in family.items()})
        _INDEX[taxid] = index
    return _INDEX[taxid]


def cai(sequence, taxid=YEAST):
    """The codon adaptation index, or None for a sequence it cannot score.

    None rather than a raise, because ytk-cds-qc measures sequences that have
    already failed their checks: a gene in the wrong frame still deserves a row
    in the table, just not a CAI on it.
    """
    if not sequence or len(sequence) % 3 or set(sequence) - set("ACGT"):
        return None
    index = _index(taxid)
    # Mirrors what calculate() counts: it skips ATG and TGG, which have no
    # synonym, and the stop codons, which are not in the index. A sequence of
    # nothing but those has no codon to score, so it has no CAI.
    if not any(sequence[at:at + 3] not in ("ATG", "TGG")
               and sequence[at:at + 3] in index
               for at in range(0, len(sequence), 3)):
        return None
    return round(index.calculate(sequence), 3)


def measure(sequence):
    """Every number, ready to be a row. Keys are the input.csv column names."""
    sequence = sequence.upper()
    low, high = gc_windows(sequence)
    base, run, _ = longest_run(sequence)
    repeat, _, _ = longest_repeat(sequence)
    return {
        "gc": gc(sequence),
        "gc_window_min": low,
        "gc_window_max": high,
        "longest_run": run,
        "longest_run_base": base,
        "longest_repeat": repeat,
        "bsmbi": enzymes.count(sequence, "BsmBI"),
        "bsai": enzymes.count(sequence, "BsaI"),
        "cai_scer": cai(sequence, YEAST),
    }


# One minimal coding sequence, only so the column order lives in one place.
COLUMNS = list(measure("ATGAAATAA"))
