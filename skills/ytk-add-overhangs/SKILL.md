---
name: ytk-add-overhangs
description: Add the MoClo Yeast Toolkit (YTK) Golden Gate flanks to a DNA sequence, so the fragment can be ordered and later cut into a YTK part plasmid. Puts the BsmBI and BsaI sites and the right four-base overhang pair around the sequence for the part type given, and writes a GenBank or SnapGene .dna file. Use this whenever someone has a sequence and wants overhangs, flanks, adapters, a sequence to order or synthesise, a gene "made YTK compatible", or asks what to send to Twist or IDT for a Type 3 coding sequence or any other YTK part type. For the whole job at once, use ytk-workflow.
---

# Add YTK overhangs to a sequence

This skill wraps a sequence in the flanks a YTK part needs. The result is the
fragment to order from a gene synthesis company.

It does not build a plasmid map. Use **ytk-clone** for that.

## What the person must tell you

**One thing is required: the sequences**, as a FASTA or CSV file. Everything
else has a default that is right nearly always.

Ask for these only if they matter to the person:

- **Where to save** the output. Without `--outdir` it goes in a `ytk_output/`
  folder beside the input file, and the script prints where.
- **A different part type.** `--type` defaults to **3**, a whole coding
  sequence, which is the only type this plugin builds. The default is printed
  on every run.
- **A different format.** Maps are GenBank, which SnapGene also opens.
  `--format dna` writes SnapGene `.dna` instead, here and for every later step:
  this is the one place the format is chosen. The synthesis order is always
  written as CSV either way.

**A part type is still never read off the DNA.** Type 3 against 3a against 3b
is a design decision. If someone says they want a fusion half, pass `--type`
yourself - do not infer it from a file name, a gene name, or the fact that a
sequence starts with ATG.

## Run it

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/add_overhangs.py" --input parts.csv
```

One sequence instead of a file: write it to a one-line CSV first, so there is
only ever one input path to reason about.

```bash
printf 'name,sequence\nPi_fim_NCS_c1,ATGGCG...TAA\n' > parts.csv
```

`--input` takes the same files as the other skills: FASTA (`.fa`, `.fasta`)
or CSV with a name column and a sequence column. One part type applies to
every row, so keep one file per part type.

## Output

Everything goes in a `fragments/` folder under the output folder:

```
fragments/
  summary.csv    name, sequence, part_type, codon_opt, notes
  <gene>.gb      one linear map per sequence
```

`sequence` is the flanked sequence: what you paste into an order form. Plasmid
names are not here. They belong with the input, and `ytk-clone` reads them from
there.

Each map is linear, and the sequence itself is labelled inside the flanks. Cut
sites are left unlabelled on purpose: SnapGene shows them live under Enzymes, so
a written-in label goes out of date. `--format dna` writes `<gene>.dna` instead.

`summary.csv` does two jobs. It is the file a synthesis order is placed from,
and it is the input to `ytk-clone`. The `name` column is the plain gene name,
because `ytk-clone` names its output files from it.

### The notes column

Per gene, what was done to it and what is still in it. **Empty when the gene was
clean and nothing was done.**

- `codon_opt (<method>)` - when the codons were optimised
- `BsmBI site removed` - per enzyme cleared, read from the `removed` column of
  the input file, which `ytk-cds-qc --remove-sites` writes
- `BsaI site in the CDS` - per enzyme still there

These are not exclusive, so a note can read
`BsmBI site removed; BsaI site in the CDS`: one enzyme can be swapped out while
another stays, because no synonymous codon removes it.

The note is built from facts each of which has its own column, and never read
back out of a sentence.

**The "still there" half is always measured again here**, on the sequence being
ordered, and nothing read from a file can suppress it. A stale column can
therefore mislabel a safe fragment, but it can never hide a site. This file
costs real money at a vendor.

The flanked sequence is always printed as well, so it can be copied straight
into SnapGene or a web form.

## The flanks

Every part gets the same outer sequence:

```
actcgacaac CGTCTCa tcGGTCTCa [left] SEQUENCE [right] tGAGACC tGAGACG gttgtggtgt
```

Reading outwards from the sequence: the BsaI site cuts the part into a later
assembly, the BsmBI site cuts it out of the ordered fragment, and the ten
outer bases are there to hold a primer.

Only the `[left]` and `[right]` bits change with the part type. They live in
`data/ytk_overhangs.tsv`, one row per type.

## Type 3 has a one-base left adapter

For type 3 the left adapter is a single `T`, not four bases. That is on
purpose. The `T` plus the gene's own `ATG` makes the four-base `TATG`
junction. So a type 3 sequence must start with `ATG`, and the script stops if
it does not.

The flank adds no stop codon. A coding sequence is expected to carry its own,
the same as in **ytk-clone**. If one does not, the script stops and writes
nothing. Ask the person whether that is intended. If it is, name that sequence
with `--no-stop-codon`:

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/add_overhangs.py" --input parts.csv --no-stop-codon Pi_fim_NCS_c1 Pi_fim_OMT_c1
```

The flag names the sequences that are exempt, one or more of them. Every other
sequence in the same run still has to carry a stop codon. The name is the one in
the `name` column of the input.

A name that is not in the input stops the run, so a typo cannot pass for an
exemption.

## Never change the sequence

Report a problem. Do not fix it.

The script stops, and writes nothing, when:

- the sequence holds anything other than A, C, G and T
- a coding sequence is not a whole number of codons
- a type needing `ATG` at the start does not have one
- the first half of a protein fusion (type 3a) ends with a stop codon
- a coding sequence has no stop codon at the end, and was not named with
  `--no-stop-codon`

It notes, and carries on, when:

- the sequence holds a BsmBI or BsaI site inside it

That goes in the `notes` column, and the same words are printed under each
sequence as it is written, so what is on screen and what gets ordered cannot
drift apart. Pass it on. An internal site does not block the ordering, and
**ytk-clone** handles it, but that part cannot be re-cut with the same enzyme
later at the bench.

Never remove an internal site with a silent codon change. That is the
person's decision, not yours.

## After this

The fragment goes to a synthesis company. When the DNA arrives, use
**ytk-clone** to build the part plasmid map.
