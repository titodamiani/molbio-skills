"""The Type IIS enzymes YTK uses, and finding their sites in a sequence.

This is the one home for the enzyme check. It used to exist twice, once in
ytk-add-overhangs and once in ytk-clone, with two different implementations.
"""

# Recognition sequence, and the same site read on the other strand. Reading
# outwards from a part: BsaI cuts the part into an assembly, BsmBI cuts the
# part out of the ordered fragment.
ENZYMES = {
    "BsmBI": ("CGTCTC", "GAGACG"),
    "BsaI": ("GGTCTC", "GAGACC"),
}

# A correctly flanked fragment carries one site per strand for each enzyme.
SITES_PER_FLANKED_FRAGMENT = 2


def count(sequence, enzyme):
    """How many times enzyme cuts, counting both strands."""
    sequence = sequence.upper()
    return sum(sequence.count(site) for site in ENZYMES[enzyme])


def positions(sequence, enzyme):
    """Where enzyme cuts, as (0-based position, the site as found)."""
    sequence = sequence.upper()
    found = []
    for site in ENZYMES[enzyme]:
        start = sequence.find(site)
        while start != -1:
            found.append((start, site))
            start = sequence.find(site, start + 1)
    return sorted(found)


def in_sequence(sequence):
    """Enzymes that cut inside a bare sequence, which has no flanks yet.

    Any site at all is a problem here: a coding sequence is meant to carry
    none, and one inside it breaks the Golden Gate assembly.
    """
    return [enzyme for enzyme in ENZYMES if count(sequence, enzyme)]


def beyond_flanks(fragment):
    """Enzymes that cut a flanked fragment somewhere other than its flanks.

    The flanks contribute two sites per enzyme by design, so only the extras
    are reported.
    """
    return [enzyme for enzyme in ENZYMES
            if count(fragment, enzyme) > SITES_PER_FLANKED_FRAGMENT]
