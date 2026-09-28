---
name: ytk-annotate-map
description: Put labels on plasmid map files - SnapGene .dna or GenBank - using the published MoClo Yeast Toolkit (YTK) parts table. Finds ColE1, CamR, the CamR promoter and terminator, any YTK part, and the cloned gene, then writes the features back into the file. Use this whenever someone wants features, labels or annotations added to a plasmid map, says a .dna file looks blank or empty in SnapGene, says the features are missing or wrong, or asks what is actually in a plasmid map.
---

# Label a plasmid map

This skill reads a `.dna` file, looks for known YTK parts in it, and writes
the file back with those parts labelled.

## Run it

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/annotate_map.py" FILE.dna
```

Several files at once, and with the cloned genes labelled too:

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/annotate_map.py" OUTDIR/*.dna --genes GENES
```

`$CLAUDE_SKILL_DIR` is this skill's own folder, so the command works from any
directory. The script finds the shared parts table by itself.

The file format is taken from the extension, and the file is written back in the
same format. A `.gb` or `.gbk` is read and written as GenBank, a `.dna` as
SnapGene. There is no flag to set.

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
inherits the error. The wrong label looks normal, so nobody notices it.

## How parts are found

Each part is found by searching for its sequence in the plasmid, treated as a
circle. That means the labels stay right no matter where the map starts, and
a part that runs past the end of the circle is still labelled correctly.

If a part sequence turns up more than once the script stops, because then any
single answer would be a guess.

## What is not labelled

Cut sites and fusion scars are left off on purpose. SnapGene shows those live
under Enzymes and Common Features. A written-in label goes out of date as soon
as the map changes, and a wrong label is worse than no label.

## After this

Use **ytk-verify-map** to check the results.
