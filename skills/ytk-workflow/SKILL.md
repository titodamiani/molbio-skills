---
name: ytk-workflow
description: Run the whole MoClo Yeast Toolkit (YTK) job for one coding sequence or a batch of them - check the sequences, add the Golden Gate flanks, design the PCR primers, clone into pYTK001, label the maps and check them. Use this whenever someone wants the complete job rather than one step, has a file of genes and wants finished primers and labelled SnapGene maps, or says something like "do the whole thing", "the full run", "all of it", "start to finish", "everything for these genes", or "just give me the primers and the maps".
---

# The whole YTK job

Runs the six skills in order, on one sequence or on a batch.

| Step | Skill | Gives |
|---|---|---|
| 1 | `ytk-cds-qc` | are these valid Type 3 coding sequences? |
| 2 | `ytk-add-overhangs` | the flanked fragment, full 10-base pad |
| 3 | `ytk-design-primers` | one CSV of primer pairs, 4-base pad |
| 4 | `ytk-clone` | the part plasmid map |
| 5 | `ytk-annotate-map` | the same maps, labelled |
| 6 | `ytk-verify-map` | the maps checked |

## Before you start

Ask for, and never guess:

- **the gene file** (`.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`)
- **where to put the output**
- **the plasmid numbers**, if they want particular ones. Real plasmid numbers
  are not in order, so they go in a `plasmid` column in the CSV. Never invent one.
- **the map format**, if they have a preference. SnapGene `.dna` is the default.
  GenBank is for anyone without SnapGene. Ask only if it might matter to them.

The part type is always 3. There is nothing to ask about that.

## Run it

**Step 1 — check the sequences.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-cds-qc/scripts/cds_qc.py" --input genes.csv
```

If anything fails here, stop and show the table. Do not carry on with a
sequence that is not a valid Type 3 CDS.

**Step 2 — add the flanks.** `--format csv` so the next step can read it.

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-add-overhangs/scripts/add_overhangs.py" \
    --input genes.csv --type 3 --format csv --outdir out/
```

**Step 3 — design the primers.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-design-primers/scripts/design_primers.py" \
    --input genes.csv --outdir out/
```

**Step 4 — clone.** Reads the flanked fragments from step 2.

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-clone/scripts/clone.py" \
    --input out/ytk_order.csv --outdir out/
```

Step 4 writes the maps to `out/maps/` and the linear fragments to
`out/fragments/`, so the next two steps can just take `out/maps/*`. Plasmid
names come from the `plasmid` column and are carried through, so the maps are
named `pTP412.dna` and not after the gene.

For GenBank, add `--format genbank` to steps 2 and 4. The maps are then
`out/maps/*.gb`, and steps 5 and 6 read them without any extra flag. Use the
same format for both steps, so the whole run is in one format.

**Step 5 — label the maps.** `--genes` points at the **original** gene file, so
the gene itself gets labelled and not the flanked version.

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-annotate-map/scripts/annotate_map.py" \
    out/maps/*.dna --genes genes.csv
```

**Step 6 — check the maps.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-verify-map/scripts/verify_map.py" \
    out/maps/*.dna --genes genes.csv
```

With no reference file this checks the map on its own: reading frame, ATG, stop
codon, and that the start of the map does not split the gene. If they have
hand-made maps to compare against, add `--reference`.

## Stop on a problem

Steps 1, 3 and 4 can stop. When one does:

- **Show the table.** Every skill collects a whole batch's problems and prints
  one table. Do not paraphrase it.
- **Ask once**, not once per sequence.
- **Never fix a sequence to get past a check.** A cut site inside a coding
  sequence has to be removed deliberately, and only when asked.

## Report at the end

- how many pairs of primers, and where the CSV is
- how many maps, and where
- which genes hold an internal BsmBI or BsaI site, because those cannot be
  re-cut with that enzyme later
- any primer pair that would not balance within 2 C, which needs finishing by
  hand
- anything skipped, and why

## The two outputs that differ on purpose

Step 2 keeps a **10-base pad**. Step 3 trims it to **4**. That is not a bug: the
pad sits outside both BsmBI sites, so it is cut off and thrown away, and both
give the same plasmid. Step 2's fragment is what you order from a synthesis
company if the PCR keeps failing.
