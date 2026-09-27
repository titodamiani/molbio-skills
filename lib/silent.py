"""Remove a cut site from a coding sequence without changing the protein.

Golden Gate cuts a BsmBI or BsaI site wherever it finds one, so a site inside a
coding sequence takes the part apart. It can be removed by swapping one codon
for a different codon of the same amino acid.

**This is the only code here that changes a sequence, and it runs only when
asked for by name.** Nothing calls it on its own. Every changed base is
reported, and the protein is translated before and after and compared: if it
differs by one residue the change is thrown away.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cds
import deps
import enzymes

deps.require("Bio")
from Bio.Seq import Seq  # noqa: E402

BASES = "ACGT"


def _synonyms():
    """Codons grouped by the amino acid they make."""
    groups = {}
    for first in BASES:
        for second in BASES:
            for third in BASES:
                codon = first + second + third
                groups.setdefault(str(Seq(codon).translate()), []).append(codon)
    return groups


SYNONYMS = _synonyms()


def protein(sequence):
    return str(Seq(sequence).translate())


def sites(sequence):
    """Every cut site, as (position, length, enzyme). Both strands."""
    found = []
    for name in enzymes.ENZYMES:
        for position, site in enzymes.positions(sequence, name):
            found.append((position, len(site), name))
    return sorted(found)


def _codons_touching(position, length):
    """The codon start positions a site overlaps."""
    first = (position // 3) * 3
    last = ((position + length - 1) // 3) * 3
    return range(first, last + 1, 3)


def _differences(before, after):
    return [index for index, (a, b) in enumerate(zip(before, after)) if a != b]


def _swap(sequence, at, codon):
    return sequence[:at] + codon + sequence[at + 3:]


def _best_change(sequence, position, length):
    """The smallest codon swap that removes this site and adds none.

    Fewest changed bases wins, then a change at the third base of the codon,
    which is where the genetic code is most often redundant.
    """
    before = len(sites(sequence))
    options = []
    for at in _codons_touching(position, length):
        codon = sequence[at:at + 3]
        if len(codon) < 3:
            continue
        amino = protein(codon)
        # A stop codon is the terminator, and Met has no synonym anyway.
        if amino in ("*", "M"):
            continue
        for candidate in SYNONYMS.get(amino, []):
            if candidate == codon:
                continue
            changed = _swap(sequence, at, candidate)
            if len(sites(changed)) >= before:
                continue
            moved = _differences(codon, candidate)
            options.append((len(moved), -max(moved), at, codon, candidate))
    if not options:
        return None
    _, _, at, codon, candidate = min(options)
    return at, codon, candidate


def remove_sites(sequence):
    """A sequence with no internal cut sites, and the list of changes.

    Returns (new_sequence, changes, remaining). `changes` is one entry per base
    moved. `remaining` lists any site that could not be removed, so a caller can
    say so rather than pretend the job is done.
    """
    sequence = sequence.upper()
    faults = cds.problems(sequence)
    if faults:
        raise ValueError("not a Type 3 coding sequence: " + "; ".join(faults))

    original = sequence
    wanted = protein(original)
    changes = []
    # Every pass removes at least one site, so this cannot run away.
    while True:
        remaining = sites(sequence)
        if not remaining:
            break
        change = _best_change(sequence, remaining[0][0], remaining[0][1])
        if change is None:
            break
        at, codon, candidate = change
        sequence = _swap(sequence, at, candidate)
        for offset in _differences(codon, candidate):
            changes.append({
                "position": at + offset + 1,
                "was": codon[offset],
                "now": candidate[offset],
                "codon": f"{codon} -> {candidate}",
                "amino_acid": protein(codon),
                "codon_number": at // 3 + 1,
            })

    if protein(sequence) != wanted:
        raise AssertionError(
            "the protein changed, so nothing was kept. This is a bug in "
            "lib/silent.py, not in your sequence")
    if len(sequence) != len(original):
        raise AssertionError("the length changed, so nothing was kept")

    return sequence, changes, sites(sequence)


def describe(changes):
    """The changes as a table, one row per base."""
    if not changes:
        return "no bases changed"
    lines = [f"{'base':>7s} {'codon':>6s} {'was':>4s} {'now':>4s}  "
             f"{'codon change':14s} amino acid"]
    for change in changes:
        lines.append(f"{change['position']:7d} {change['codon_number']:6d} "
                     f"{change['was']:>4s} {change['now']:>4s}  "
                     f"{change['codon']:14s} {change['amino_acid']}")
    return "\n".join(lines)
