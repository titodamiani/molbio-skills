---
name: ytk-verify-map
description: Check SnapGene .dna files. Compares a map against a reference as a circle, so a different start point is not reported as a difference, and checks the reading frame, the linear or circular flag, the double-strand flag, and that the start of the map does not split the gene. Use this whenever someone wants to verify, compare or sanity-check plasmid maps or .dna files, asks whether a generated map matches one made by hand, asks why two maps of the same plasmid look different, or says a file opens oddly in SnapGene.
---

# Check plasmid files

This skill checks `.dna` files and prints one line per check.

## Run it

Check a file on its own:

```bash
python3 scripts/verify_map.py FILE.dna --genes GENES
```

Compare one file against a reference:

```bash
python3 scripts/verify_map.py FILE.dna --reference REFERENCE.dna --genes GENES
```

`scripts/verify_map.py` sits in this skill's own folder. `--reference` takes one file
at a time. Exit code 1 means something failed.

## Compare circles, not text

A plasmid is a loop. Two files can hold exactly the same DNA while starting at
different points on that loop. Comparing them as plain text would report a
difference that is not there.

So the check turns one sequence until it lines up with the other. The output
says either `matches the reference exactly` or `matches the reference, turned
by N bp`. Both mean the DNA is the same.

This comes up often with files made by hand in SnapGene, because a person
picks whatever start point the software offered that day.

## What gets checked

Every file:

- the double-strand bit is set — with it off, SnapGene treats the file as a
  single strand and hides the Enzymes tab
- the length

Circular files, with `--genes`:

- the gene is not split by the start of the map
- the gene is a whole number of codons
- the gene starts with ATG
- the gene ends with a stop codon

Linear files, with `--genes`:

- the sequence is unchanged from the input

With `--reference`, also whether the two describe the same DNA.

## Reading the result

A failure is a real problem and is worth stopping for. Report exactly which
check failed and for which file. Do not paper over it, and do not edit
anyone's sequence to make a check pass — the check is there to catch the
mistake, not to be satisfied.

A `turned by N bp` line is not a failure. Say so plainly if someone asks,
because it reads like one.
