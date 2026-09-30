"""The notes column: one sentence saying what was done to a sequence, and what
is still in it.

    codon_opt (JCat); BsmBI site removed; BsaI site in the CDS

Built from facts, and never read back out of a cell. Every file that carries a
note carries the facts behind it in their own columns - `removed` and
`codon_opt_method` - so a later step rebuilds the sentence instead of parsing
it. That is why there is no sentinel here for "nothing to say": an empty note is
an empty string.

The clauses are not exclusive. A gene can have its BsmBI site swapped out while
a BsaI site stays put, because no synonymous codon removes it, so a note joins
them rather than picking one.
"""
# Inside the `removed` column, which holds enzyme names and not prose.
SEPARATOR = ";"


def summarise(removed, present, codon_opt_method=None):
    """The notes cell for one sequence.

    `removed` and `present` are collections of enzyme names. `removed` has to
    travel in its own column, because a site that was swapped out cannot be
    measured again afterwards. `present` is always measured fresh on the
    sequence in hand, so a stale column can mislabel a safe fragment but can
    never hide a site and bless an unsafe one.
    """
    clauses = []
    if codon_opt_method:
        clauses.append(f"codon_opt ({codon_opt_method})")
    clauses += [f"{enzyme} site removed" for enzyme in sorted(removed)]
    clauses += [f"{enzyme} site in the CDS" for enzyme in sorted(present)]
    return "; ".join(clauses)


def pack(enzyme_names):
    """Enzyme names as the `removed` column holds them."""
    return SEPARATOR.join(sorted(enzyme_names))


def unpack(cell):
    """The `removed` column back as a set of enzyme names."""
    return {name.strip() for name in cell.split(SEPARATOR) if name.strip()}
