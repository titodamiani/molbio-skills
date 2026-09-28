---
name: ytk-add-overhangs
description: Add the MoClo Yeast Toolkit (YTK) Golden Gate flanks to a DNA sequence, so the fragment can be ordered and later cut into a YTK part plasmid. Puts the BsmBI and BsaI sites and the right four-base overhang pair around the sequence for the part type given, and writes a SnapGene .dna or GenBank file. Use this whenever someone has a sequence and wants overhangs, flanks, adapters, a sequence to order or synthesise, a gene "made YTK compatible", or asks what to send to Twist or IDT for a Type 3 coding sequence or any other YTK part type. For the whole job at once, use ytk-workflow.
---

# Add YTK overhangs to a sequence

This skill wraps a sequence in the flanks a YTK part needs. The result is the
fragment to order from a gene synthesis company.

It does not build a plasmid map. Use **ytk-clone** for that.

## What the person must tell you

Three things are required. Never guess any of them.

1. **The sequence.** Plain DNA, or a FASTA or CSV file.
2. **The part type.** For example type 3 for a coding sequence. The overhang
   pair depends only on this. The sequence does not say what type it is, so a
   guess here silently puts the part in the wrong slot of an assembly.
3. **The name** of the sequence, for example `Pi_fim_NCS_c1`. The output file
   is named after it.

Ask for these two as well, if they matter to the person:

4. **Where to save** the output. The default is the folder you are in.
5. **A different format.** The default output is a SnapGene `.dna` file.
   GenBank is there for anyone without SnapGene. CSV and FASTA are there for
   sending sequences to a synthesis company.

If the part type is missing, ask. Do not read it off a file name, a gene name
or the fact that a sequence starts with ATG.

## Run it

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/add_overhangs.py" --sequence ATGGCG... --type 3 --name Pi_fim_NCS_c1
```

For a batch:

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/add_overhangs.py" --input parts.csv --type 3 --outdir out/
```

`--input` takes the same files as the other skills: FASTA (`.fa`, `.fasta`)
or CSV with a name column and a sequence column. One `--type` applies to
every row, so keep one file per part type.

## Output

One map per sequence, called `<name>_oh.dna`. It is linear, and the sequence
itself is labelled inside the flanks. Cut sites are left unlabelled on purpose.
SnapGene shows them live under Enzymes, so a written-in label goes out of date.

`--format genbank` writes `<name>_oh.gb` instead of `.dna`.
`--format fasta` writes `<name>_oh.fasta` instead, one file per sequence.
`--format csv` writes a single `ytk_order.csv` for the whole batch, with one
row per part. Use CSV for an order form.

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
python3 "$CLAUDE_SKILL_DIR/scripts/add_overhangs.py" --input parts.csv --type 3 --no-stop-codon Pi_fim_NCS_c1 Pi_fim_OMT_c1
```

The flag names the sequences that are exempt, one or more of them. Every other
sequence in the same run still has to carry a stop codon. For a single
sequence the name is whatever `--name` says, so it reads
`--name demo --no-stop-codon demo`.

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

It warns, and carries on, when:

- the sequence holds a BsmBI or BsaI site inside it

Pass that on. An internal site does not block the ordering, and **ytk-clone**
handles it, but that part cannot be re-cut with the same enzyme later at the
bench.

Never remove an internal site with a silent codon change. That is the
person's decision, not yours.

## Part types we are not sure about

Some rows of `data/ytk_overhangs.tsv` carry a note in the `unsure` column. For
those part types the script prints the note, stops, and writes nothing.

Show the message to the person and ask whether to go ahead. Only if they say
yes, run it again with `--approve-unsure`. Never add that flag on your own
initiative — the whole point of the gate is that a person decides.

When the answer is confirmed, clearing that cell in the table turns the gate
off. No code changes.

## After this

The fragment goes to a synthesis company. When the DNA arrives, use
**ytk-clone** to build the part plasmid map.
