---
name: ytk-annotate
description: Put labels on SnapGene .dna plasmid maps using the published MoClo Yeast Toolkit (YTK) parts table. Finds ColE1, CamR, the CamR promoter and terminator, any YTK part, and the cloned gene, then writes the features back into the file. Use this whenever someone wants features, labels or annotations added to a plasmid map, says a .dna file looks blank or empty in SnapGene, says the features are missing or wrong, or asks what is actually in a plasmid map.
---

# Label a plasmid map

This skill reads a `.dna` file, looks for known YTK parts in it, and writes
the file back with those parts labelled.

## Run it

```bash
python3 scripts/annotate.py FILE.dna
```

Several files at once, and with the cloned genes labelled too:

```bash
python3 scripts/annotate.py OUTDIR/*.dna --genes GENES
```

`scripts/annotate.py` sits in this skill's own folder. It finds the shared
parts table by itself, so it works wherever the plugin is installed.

`--genes` takes the same FASTA or CSV file used to build the plasmids. Each
gene in that file is searched for, and the one that is present gets its own
name on the map. No pairing of gene to plasmid is needed.

## Where labels come from

Labels come from `data/ytk_parts.tsv` at the root of the plugin. That table
was built from the published GenBank files of the toolkit (Lee et al., ACS
Synthetic Biology 2015).

Never copy labels from another map file, even one that looks right. Map files
carry mistakes. One known case: a `pYTK001.dna` in circulation says the ColE1
origin covers the whole plasmid. It is 764 bp. Anything built on that file
inherits the error, and it is hard to spot later because the label looks
perfectly ordinary.

## How parts are found

Each part is found by searching for its sequence in the plasmid, treated as a
circle. That means the labels stay right no matter where the map starts, and
a part that runs past the end of the circle is still labelled correctly.

If a part sequence turns up more than once the script stops, because then any
single answer would be a guess.

## What is not labelled

Cut sites and fusion scars are left off on purpose. SnapGene shows those live
under Enzymes and Common Features. A fixed label for them would go stale as
soon as anything changed, and a stale label is worse than none.

## After this

Use the **ytk-qc** skill to check the results.
