"""Rewrite a coding sequence with the most-used codon of every amino acid.

The second piece of code here that changes a sequence, after lib/silent.py, and
like that one it runs only when asked for by name: only ytk-codon-optimise calls
it. Nothing calls it on its own. The protein is translated before and after and
compared, and if it differs at all nothing is kept.

The usage figures come from python_codon_tables, which gives each amino acid's
codons as a share of that amino acid's own total, so the most-used codon is the
largest share and no arithmetic is needed to find it.

Biopython has nearly this in Bio.SeqUtils.CodonAdaptationIndex.optimize(), and
it is not used for two reasons: it swaps the stop codon for whichever stop is
most used, where a gene here keeps its own, and on a tie it either raises or
picks silently rather than by a stated rule. tests/test_codons.py pins this
module against it on the codons where the two should agree, so the borrowed
logic cannot drift unnoticed.

**Nothing here avoids BsmBI or BsaI sites, and a rewrite can create one.**
Removing them is ytk-remove-cut-sites' job, and it has to run after this.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cds
import deps
import silent

deps.require("python_codon_tables")
import python_codon_tables as pct  # noqa: E402

# Saccharomyces cerevisiae. Any taxid the package carries works.
YEAST = 4932

# How many times in a row the same codon may be used. A synthesis company can
# refuse an order over a repeat, and the argmax rule builds one whenever the
# protein repeats an amino acid: poly-Q becomes CAACAACAACAA. Capping the run at
# two and using the runner-up for the third leaves every doubled residue alone,
# which is most of them, and only touches a genuine low-complexity tract.
MAX_SAME_CODON = 2

_CHOICES = {}


def ranked(taxid=YEAST):
    """{amino acid: its codons, most used first}.

    Ties are broken alphabetically by codon. The shares in the table are
    rounded to two places, so two codons really can tie, and without a stated
    rule the winner would follow whatever order the table happened to be built
    in - which makes the same gene recode differently on another machine.
    """
    table = pct.get_codons_table(taxid)
    return {amino: sorted(shares, key=lambda codon: (-shares[codon], codon))
            for amino, shares in table.items()}


def method(taxid=YEAST):
    """What went in the codon_opt_method column.

    The package version is part of it because the table lives in the package: a
    bump that moves one codon would otherwise change a sequence with nothing in
    the output saying why.
    """
    return f"argmax/{taxid}/python_codon_tables-{pct.__version__}"


def _choices(taxid):
    """{amino acid: (most used codon, the codon to break a run with)}.

    The second entry is the runner-up, or the same codon again for Met and Trp,
    which have no second. Which one gets used is decided in recode, by how many
    of the first are already there.
    """
    if taxid not in _CHOICES:
        _CHOICES[taxid] = {
            amino: (order[0], order[1] if len(order) > 1 else order[0])
            for amino, order in ranked(taxid).items()}
    return _CHOICES[taxid]


def recode(sequence, taxid=YEAST):
    """A sequence of most-used codons, and the list of changes.

    One change per codon moved, not per base. lib/silent.py reports per base
    because it moves one; a rewrite moves most of the gene, and a per-base list
    of that is noise rather than a record.
    """
    sequence = sequence.upper()
    # "full": swapping a codon needs a reading frame, so this only ever works on
    # a whole coding sequence.
    faults = cds.problems(sequence, "full")
    if faults:
        raise ValueError("not a Type 3 coding sequence: " + "; ".join(faults))
    # A sequence read in the wrong frame can still be a whole number of codons
    # that starts ATG and ends in a stop, and it rewrites into something that
    # looks clean and means nothing. An internal stop is the one cheap sign of
    # that, so here it stops the run, where ytk-cds-qc only warns about it.
    position = cds.internal_stop(sequence)
    if position is not None:
        raise ValueError(f"has a stop codon at base {position + 1}, before the "
                         "end, so the reading frame is probably wrong")

    choices = _choices(taxid)
    wanted = silent.protein(sequence)
    # The first codon is the ATG and the last is the gene's own stop. Neither is
    # a choice, so the loop covers only what lies between them.
    out = [sequence[:3]]
    changes = []
    for at in range(3, len(sequence) - 3, 3):
        codon = sequence[at:at + 3]
        amino = silent.protein(codon)
        best, runner_up = choices[amino]
        run = out[-MAX_SAME_CODON:]
        repeated = len(run) == MAX_SAME_CODON and set(run) == {best}
        out.append(runner_up if repeated else best)
        if out[-1] != codon:
            changes.append({
                "codon_number": at // 3 + 1,
                "amino_acid": amino,
                "was": codon,
                "now": out[-1],
            })
    out.append(sequence[-3:])
    new = "".join(out)

    if silent.protein(new) != wanted:
        raise AssertionError(
            "the protein changed, so nothing was kept. This is a bug in "
            "lib/codons.py, not in your sequence")
    if len(new) != len(sequence):
        raise AssertionError("the length changed, so nothing was kept")

    return new, changes


def describe(changes):
    """The changes as a table, one row per codon."""
    if not changes:
        return "no codons changed"
    lines = [f"{'codon':>6s} {'was':>4s} {'now':>4s}  amino acid"]
    for change in changes:
        lines.append(f"{change['codon_number']:6d} {change['was']:>4s} "
                     f"{change['now']:>4s}  {change['amino_acid']}")
    return "\n".join(lines)
