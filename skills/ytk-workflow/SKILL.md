---
name: ytk-workflow
description: Run the whole MoClo Yeast Toolkit (YTK) job for one coding sequence or a batch of them - check the sequences, add the Golden Gate flanks, design the PCR primers, clone into an entry vector, label the maps and check them. Use this whenever someone wants the complete job rather than one step, has a file of genes and wants finished primers and labelled plasmid maps, or says something like "do the whole thing", "the full run", "all of it", "start to finish", "everything for these genes", or "just give me the primers and the maps".
---

# The whole YTK job

Runs the six skills in order, on one sequence or on a batch.

| Step | Skill | Gives |
|---|---|---|
| 1 | `ytk-cds-qc` | are these valid Type 3 coding sequences? |
| 2 | `ytk-add-overhangs` | the flanked fragments, and the synthesis order |
| 3 | `ytk-design-primers` | one CSV of primer pairs, 4-base pad |
| 4 | `ytk-clone` | the part plasmid maps |
| 5 | `ytk-annotate-map` | the same maps, labelled |
| 6 | `ytk-verify-map` | the maps checked |

## Before you start

Ask for, and never guess:

- **the gene file** (`.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`)
- **where to put the output**
- **whether to remove internal BsmBI/BsaI sites.** This is the one question that
  changes the DNA. Never assume it.

Everything else has a default that is right nearly always, so do not ask about
it. The part type defaults to 3 and every step prints the type it used. Ask about
it only if the person mentions a fusion half, a promoter, a terminator, a tag or
a non-YTK vector. **Never read the part type off the DNA** - 3 against 3a against
3b is their design decision. Maps are GenBank unless they ask for SnapGene
`.dna`, and that is chosen once, at step 2. Plasmid names come from a `plasmid`
column, with `<gene>_<backbone>` as the fallback.

**A part type other than 3 skips step 1.** `ytk-cds-qc` asks one question, *is
this a whole coding sequence*, so it would refuse a promoter or a fusion half
with "not a Type 3 CDS". For any other type, start at step 2, and take the
internal-cut-site information from the `notes` column that step 2 writes.

**If they do not say where to put the output**, do not ask either. Leave
`--outdir` off and every step writes to a `ytk_output/` folder beside the input
file, printing where it went. Then tell them the path in your report.

**If they paste sequences into the chat instead of giving a file**, write them
to a CSV first and treat that as the input. Everything below then works
unchanged:

```
printf 'name,sequence\nmy_gene,ATGAAA...TAA\n' > genes.csv
```

## The prompt to hand out

When someone asks for "the prompt", "the template" or "what do I paste", print
this and nothing else:

```
Use the ytk-workflow skill to run the full YTK workflow on these genes.

Input:   /path/to/genes.csv
Output:  /path/to/output/
Part type (default 3):
Remove internal BsmBI/BsaI sites in the CDS: [True/False]
```

## What the output looks like

Always this, whatever the input was:

```
OUTPUT/
  input.csv              a copy of the input, as CSV whatever it was
  corrected_genes/       only when sites were removed
    input_corrected.csv  what every later step reads
    summary.csv          the old sequence beside the new one
  ytk_primers.csv        two rows per gene, one per oligo
  fragments/
    summary.csv          the synthesis order
    <gene>.gb            the flanked sequence, one per gene
  plasmids/
    <backbone>/
      summary.csv
      <plasmid>.gb
```

## Run it

Set these first. Do not copy the input in by hand: step 1 writes `input.csv`
itself, which turns a FASTA or a `.dna` file into the CSV every later step reads.

```
OUT=/path/to/output
TYPE=3
BACKBONE=pYTK001
ENZYME=BsmBI
MAPS="$OUT/plasmids/$BACKBONE"
```

**Step 1 — check the sequences.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-cds-qc/scripts/cds_qc.py" --input genes.csv
```

If anything fails here, stop and show the table. Do not carry on with a
sequence that is not a valid Type 3 CDS.

An internal BsmBI or BsaI site is a warning, not a failure. Show the table and
ask once whether to remove them. **If they say yes:**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-cds-qc/scripts/cds_qc.py" \
    --input genes.csv --remove-sites --outdir "$OUT"
```

That writes `$OUT/corrected_genes/`. **Every later step then reads
`input_corrected.csv` from there, not `input.csv`.** Steps 3, 5 and 6 included:
primers designed from the original sequence would not match the corrected gene,
and a map verified against the old sequence would report a difference that is
not there.

Below, `GENES` means `$OUT/input.csv` normally and
`$OUT/corrected_genes/input_corrected.csv` when the removal step ran.

**Step 2 — add the flanks.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-add-overhangs/scripts/add_overhangs.py" \
    --input "$GENES" --type "$TYPE" --outdir "$OUT"
```

Writes `$OUT/fragments/`: one `<gene>.gb` per gene, plus `summary.csv`. That
summary is both the synthesis order and the input to step 4. It picks up the
`removed` column from `$GENES` when there is one, so a gene that had a site
removed says so on the row its sequence is on.

**This is where the map format is chosen.** Add `--format dna` here for SnapGene
files, and step 4 follows by itself.

**Step 3 — design the primers.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-design-primers/scripts/design_primers.py" \
    --input "$GENES" --type "$TYPE" --outdir "$OUT"
```

**Step 4 — clone.** Reads the fragment summary from step 2, and the plasmid
names from `$GENES`.

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-clone/scripts/clone.py" \
    --input "$OUT/fragments/summary.csv" --genes "$GENES" --outdir "$OUT" \
    --type "$TYPE" --backbone "$BACKBONE" --enzyme "$ENZYME"
```

Writes `$MAPS`, named after the backbone actually used. Plasmid names come from
the `plasmid` column of `$GENES`, so a map is called `pTP412.gb` and not after
the gene.

The fragment file from step 2 carries the part type, so a `--type` that disagrees
with it stops the run before anything is written.

Do not pass `--format` here. It matches the fragments from step 2 by itself.

**Step 5 — label the maps.** `--genes` points at `$GENES`, so the gene itself
gets labelled and not the flanked version.

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-annotate-map/scripts/annotate_map.py" \
    $(find "$MAPS" -type f ! -name summary.csv) --genes "$GENES"
```

**Step 6 — check the maps.**

```
python3 "$CLAUDE_PLUGIN_ROOT/skills/ytk-verify-map/scripts/verify_map.py" \
    $(find "$MAPS" -type f ! -name summary.csv) --genes "$GENES"
```

**Do not glob `*.gb` in steps 5 and 6.** Step 4 follows step 2's format, so the
maps may be `.dna`. `summary.csv` is the only other file in that folder, so leave
it out rather than guessing the suffix.

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
- which genes had a site removed, and which still hold one, straight from the
  `notes` column of `fragments/summary.csv`. An empty note means the gene was
  clean and nothing was done to it.
- any primer pair that would not balance within 2 C, which needs finishing by
  hand
- which part type, backbone and enzyme were used
- anything skipped, and why

## The two outputs that differ on purpose

Step 2 keeps a **10-base pad**. Step 3 trims it to **4**. That is not a bug: the
pad sits outside both BsmBI sites, so it is cut off and thrown away, and both
give the same plasmid. Step 2's fragment is what you order from a synthesis
company if the PCR keeps failing.
