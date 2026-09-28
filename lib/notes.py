"""The notes column: what was cleared out of a CDS, and what is still in it.

One module, because the wording travels between two skills. `ytk-cds-qc` writes
it when it removes a cut site, and `ytk-add-overhangs` carries it into the order
file. If the two ever disagreed about the phrasing, the merge in `for_order`
would silently drop half of it.

The two facts are **not** exclusive. A gene can have its BsmBI site swapped out
while a BsaI site stays put, because no synonymous codon removes it. So a note
joins them rather than picking one:

    BsmBI site removed; BsaI site in the CDS
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import enzymes

NO_SITES = "input CDS contained no inner BsmBI/BsaI cut sites"
REMOVED = "site removed"
PRESENT = "site in the CDS"


def _join(clauses):
    return "; ".join(clauses) if clauses else NO_SITES


def after_removal(before, remaining):
    """The note for a sequence that has just been through silent.remove_sites.

    `before` and `remaining` are sets of enzyme names. What was cleared is what
    was there minus what is left, because remove_sites reports which bases moved
    and which sites survive, never which enzyme it cleared.
    """
    return _join([f"{enzyme} {REMOVED}" for enzyme in sorted(before - remaining)]
                 + [f"{enzyme} {PRESENT}" for enzyme in sorted(remaining)])


def for_order(sequence, inherited):
    """The note for the order file, from the input's note plus the sequence.

    What is still in the CDS is never inherited: it is measured again, here, on
    the sequence actually being ordered. A stale or hand-edited notes column can
    therefore mislabel a safe fragment, but it can never hide a site and bless an
    unsafe one. This file is what a synthesis order is placed from, so that
    asymmetry is the whole point.

    Everything else in the cell is carried through word for word, because the
    notes column is a CSV cell someone can edit in Excel. Matching it against a
    fixed phrase would drop `BsmBI site removed.` for its full stop, and the
    all-clear would then be printed over a gene that really was changed. Only an
    empty cell can produce the all-clear.
    """
    kept = []
    for clause in inherited.split(";"):
        clause = clause.strip()
        if not clause or clause == NO_SITES:
            continue
        if clause.endswith(PRESENT):
            continue
        kept.append(clause)
    kept += [f"{enzyme} {PRESENT}" for enzyme in enzymes.in_sequence(sequence)]
    return _join(kept)
