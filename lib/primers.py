"""Design PCR primers that carry the YTK flanks.

The primer amplifies a coding sequence from cDNA and comes out already
YTK-compatible, so the product drops straight into the entry vector.

    forward:  pad CGTCTC atc GGTCTC a T    binding region
    reverse:  pad CGTCTC a   GGTCTC a GGAT binding region

Only the binding region sticks to the template in the first PCR cycle, so every
number reported here - Tm, GC, length - is measured on the binding region alone.
"""
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import deps
import enzymes
import flanks

deps.require("primer3")
import primer3  # noqa: E402

# Phusion's reaction conditions: 50 mM monovalent, 1.5 mM Mg2+, 0.2 mM dNTPs,
# 500 nM primer, nearest-neighbour thermodynamics from SantaLucia.
TM_ARGS = dict(mv_conc=50, dv_conc=1.5, dntp_conc=0.2, dna_conc=500,
               tm_method="santalucia", salt_corrections_method="schildkraut")

# Measured on 2026-09-27 against the "Tm Phusion" column of oligo_stocks, which
# was filled in from the NEB Phusion calculator by hand. Fitted over the 45 of
# 49 rows whose own stated binding length matches their sequence; TD014R,
# TD020R, TD023F and TD031R disagree with themselves and were left out. After
# this offset every one of the 45 lands within 0.51 C of the NEB value.
TM_OFFSET = 1.31

MIN_BINDING = 18          # we amplify from cDNA, so never shorter
TARGET_BINDING = (20, 25)
TARGET_LENGTH = 50        # a target, not a wall: longer just costs more
TM_FLOOR = 50             # hard
TM_LOW = 55               # below this, warn
TM_TARGET = (60, 68)
GC_TARGET = (40.0, 60.0)
GC_LOW = 30.0             # below this, warn
MAX_PAIR_GAP = 2.0
LONGEST_RUN = 4

# NEB's rule for primers over 20 nt: anneal 3 C above the lower primer Tm.
ANNEALING_MARGIN = 3


def melting_temp(binding):
    """Tm of the binding region, matched to the NEB Phusion calculator."""
    return primer3.calc_tm(binding, **TM_ARGS) + TM_OFFSET


def gc_percent(binding):
    return 100.0 * sum(binding.count(base) for base in "GC") / len(binding)


def annealing_temp(forward_tm, reverse_tm):
    """The temperature to anneal the pair at.

    Not a melting temperature, and not an average: it sits above both primer
    Tms. This reproduces the "Combined Tm Phusion" column exactly, in all 12
    rows that have one.
    """
    return round(min(forward_tm, reverse_tm)) + ANNEALING_MARGIN


def longest_run(binding):
    return max((len(run) for run in re.findall(r"(A+|C+|G+|T+)", binding)),
               default=0)


def _candidates(template, constant_length, pad):
    """Usable binding regions, as (sequence, Tm) pairs.

    A candidate must be long enough to be specific on cDNA, end in G or C, hold
    its Tm above the floor, and leave the primer inside the length target.
    """
    longest = min(TARGET_LENGTH - constant_length - pad, len(template))
    usable = []
    for n in range(MIN_BINDING, longest + 1):
        if template[n - 1] not in "GC":
            continue
        binding = template[:n]
        tm = melting_temp(binding)
        if tm >= TM_FLOOR:
            usable.append((binding, tm))
    return usable


def _distance_from_window(tm):
    """How far a Tm sits outside the target window. Zero when inside."""
    return max(0.0, TM_TARGET[0] - tm, tm - TM_TARGET[1])


def _rank(forward, reverse):
    """How good a pair is, smaller being better.

    Getting both Tms into the window comes first, because a primer that anneals
    badly does not work at all. Then the shorter pair wins, because a shorter
    oligo costs less. The gap only breaks remaining ties, since every pair that
    reaches this point already sits inside the limit.

    Distances and gaps are rounded to half a degree. Without that, a hundredth
    of a degree would decide and the length preference would never be reached.
    """
    (_, tm_f), (_, tm_r) = forward, reverse
    outside = _distance_from_window(tm_f) + _distance_from_window(tm_r)
    return (_half(outside), len(forward[0]) + len(reverse[0]),
            _half(abs(tm_f - tm_r)))


def _half(value):
    return round(value * 2) / 2


def design(name, sequence, adapters, pad=flanks.MIN_PAD):
    """One primer pair for one coding sequence.

    Returns a dict on success, or a dict with a "blocked" key when no pair can
    be made under the rules. Nothing about the sequence is ever changed.
    """
    sequence = sequence.upper()
    reverse_template = flanks.reverse_complement(sequence)
    forward_constant, reverse_constant = flanks.constant_length(adapters, pad)

    forwards = _candidates(sequence, forward_constant - pad, pad)
    reverses = _candidates(reverse_template, reverse_constant - pad, pad)
    if not forwards or not reverses:
        side = "forward" if not forwards else "reverse"
        return {"name": name,
                "blocked": f"no {side} binding region of at least {MIN_BINDING} bp "
                           f"ends in G or C and reaches {TM_FLOOR} C within "
                           f"{TARGET_LENGTH} bp"}

    # The pair must be balanced to within MAX_PAIR_GAP, so that is a filter and
    # not a preference. When nothing fits, say how close the best pair got and
    # let the person decide, rather than quietly handing back a pair that breaks
    # the rule.
    pairs = [(f, r) for f in forwards for r in reverses]
    balanced = [pair for pair in pairs if abs(pair[0][1] - pair[1][1]) <= MAX_PAIR_GAP]
    if not balanced:
        closest = min(pairs, key=lambda pair: abs(pair[0][1] - pair[1][1]))
        gap = abs(closest[0][1] - closest[1][1])
        return {"name": name,
                "blocked": f"the closest pair is {gap:.1f} C apart, over the "
                           f"{MAX_PAIR_GAP:.0f} C limit",
                "closest": _describe(name, sequence, closest, adapters, pad)}

    best = min(balanced, key=lambda pair: _rank(*pair))
    return _describe(name, sequence, best, adapters, pad)


def _describe(name, sequence, pair, adapters, pad):
    (forward_binding, tm_f), (reverse_binding, tm_r) = pair
    forward = flanks.forward_primer(forward_binding, adapters, pad)
    reverse = flanks.reverse_primer(reverse_binding, adapters, pad)

    warnings = []
    for enzyme in enzymes.in_sequence(sequence):
        warnings.append(f"the sequence holds a {enzyme} site")
    for label, binding, primer, tm in (("forward", forward_binding, forward, tm_f),
                                       ("reverse", reverse_binding, reverse, tm_r)):
        gc = gc_percent(binding)
        if tm < TM_LOW:
            warnings.append(f"{label} Tm is {tm:.0f} C, below {TM_LOW}")
        # Between GC_LOW and the top of the ideal band there is nothing to say:
        # 30-40% is what a real coding sequence usually gives and it works.
        if gc < GC_LOW:
            warnings.append(f"{label} GC is {gc:.0f}%, below {GC_LOW:.0f}")
        elif gc > GC_TARGET[1]:
            warnings.append(f"{label} GC is {gc:.0f}%, above "
                            f"{GC_TARGET[1]:.0f}%")
        if longest_run(binding) >= LONGEST_RUN:
            warnings.append(f"{label} has a run of {longest_run(binding)} "
                            "identical bases")
        if len(primer) > TARGET_LENGTH:
            warnings.append(f"{label} primer is {len(primer)} bp, over "
                            f"{TARGET_LENGTH}")
    gap = abs(tm_f - tm_r)
    if gap > MAX_PAIR_GAP:
        warnings.append(f"the pair is {gap:.1f} C apart, over {MAX_PAIR_GAP:.0f}")

    return {
        "name": name,
        "forward": forward,
        "reverse": reverse,
        "forward_binding": forward_binding,
        "reverse_binding": reverse_binding,
        "forward_tm": tm_f,
        "reverse_tm": tm_r,
        "annealing_temp": annealing_temp(tm_f, tm_r),
        "pad": pad,
        "warnings": warnings,
    }


def rows_for_csv(designed):
    """Two CSV rows per pair, matching the columns of the oligo stock sheet."""
    rows = []
    for side in ("forward", "reverse"):
        binding = designed[f"{side}_binding"]
        rows.append({
            "name": f"{designed['name']}_{side}",
            "sequence": designed[side],
            "Full length (bp)": len(designed[side]),
            "Binding region length (bp)": len(binding),
            "GC%": round(gc_percent(binding)),
            "Tm Phusion (C)": round(designed[f"{side}_tm"]),
            "Combined Tm Phusion (C)": designed["annealing_temp"],
            "warnings": "; ".join(designed["warnings"]),
        })
    return rows


CSV_COLUMNS = ["name", "sequence", "Full length (bp)",
               "Binding region length (bp)", "GC%", "Tm Phusion (C)",
               "Combined Tm Phusion (C)", "warnings"]
