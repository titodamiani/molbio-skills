---
name: ytk-cds-qc
description: Check and measure DNA sequences - one coding sequence or a batch of them - before any MoClo Yeast Toolkit (YTK) cloning is designed. Confirms each one is a valid Type 3 CDS - starts with ATG, ends with a stop codon, a whole number of codons - finds any BsmBI or BsaI site inside that would break Golden Gate, and measures GC, homopolymer runs, exact repeats and the yeast codon adaptation index into input.csv. Nothing here can change a sequence. Use this whenever someone wants a sequence checked, screened, measured or validated before cloning, asks whether a gene suits YTK or Golden Gate, asks about internal cut sites, GC content or codon usage, or says something like "are these genes ok", "check these sequences first", or "will this clone". To remove a cut site use ytk-remove-cut-sites; to rewrite the codons use ytk-codon-optimise. For a finished plasmid map use ytk-verify-map. For the whole job at once, use ytk-workflow.
---

# Check and measure a coding sequence

Reads sequences and answers two questions:

1. **Is it a valid Type 3 CDS?** Starts with `ATG`, ends with a stop codon, and
   a whole number of codons.
2. **Does it hide a cut site?** A `BsmBI` or `BsaI` site inside the sequence
   gets cut during Golden Gate, so the part comes apart.

Then it measures everything else worth knowing, and writes the numbers into
`input.csv` beside each sequence, so no later step has to measure again.

This checks **sequences**. To check a finished `.dna` plasmid map instead, use
`ytk-verify-map`.

## Run it

```
python3 "$CLAUDE_SKILL_DIR/scripts/cds_qc.py" --input genes.csv
```

One sequence instead of a file: write it to a one-line CSV first, so there is
only ever one input path to reason about.

```
printf 'name,sequence\nmy_gene,ATGAAA...TAA\n' > genes.csv
```

Input can be `.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`. For a map file,
markers and origins are set aside before the insert is picked; if more than one
feature is left, the run stops and lists them, and `--feature NAME` says which.

`input.csv` is written. **Nothing is changed, and nothing here can change a
sequence.** Exits 1 if anything failed, so it can gate a longer run.

## What it reports

A table, one row per sequence, then one table of everything wrong.

| Verdict | Meaning |
|---|---|
| `ok` | a clean Type 3 CDS with no internal cut site |
| `usable, with warnings` | a valid CDS, but something needs saying |
| `not a Type 3 CDS` | cannot be a Type 3 part as it stands |

All the faults are listed, not just the first, so a sequence can be fixed in one
go.

## The numbers, and which skill each one points at

Every one of these goes into `input.csv` on the sequence's own row, so a later
step reads the fact rather than measuring again or reading a sentence apart.

| Column | What it says | What to do |
|---|---|---|
| `bsmbi`, `bsai` | sites inside the CDS | above zero, `ytk-remove-cut-sites` |
| `cai_scer` | how close the codons are to yeast's favourites, 0 to 1 | low, and `ytk-codon-optimise` is worth offering |
| `gc` | percent G plus C | nothing here fixes it |
| `gc_window_min`, `gc_window_max` | the worst 50-base stretches | a gene can average 45% and still hold a window at 20%, which is the part a synthesis company struggles with |
| `longest_run`, `longest_run_base` | the longest run of one base | a long run can get an order refused |
| `longest_repeat` | the longest exact repeat, 15 bases or more | same |

**There is one CAI column and it is against yeast.** CAI cannot say which
organism a gene was optimised for, however much it looks like it should: it is a
geometric mean of each codon's share of its own family, so a table with more even
codon usage scores every sequence higher. A random sequence that nothing
optimised scores higher against the human table than the yeast one. Never claim a
gene "was optimised for" anything from a CAI number.

## What to do about a cut site

**Report it and stop. Do not remove it unless the person asks in that message.**

This skill has no way to remove one. `ytk-remove-cut-sites` does, and it moves
one base per site on a real gene.

## Never change a sequence to make a check pass

If a sequence fails, report it. Removing a cut site or rewriting the codons is
sometimes the right call in the lab. It is never the right call to silence a
warning.

## A stop codon in the middle

Reported as a warning, not a failure. It usually means one of two things: the
reading frame is wrong, or the sequence is not the coding sequence it was
thought to be. Say this before anyone orders DNA.

It is a warning here and a refusal in `ytk-codon-optimise`, because a sequence in
the wrong frame rewrites into something that looks clean and means nothing.

## After this

- `ytk-remove-cut-sites` if a site has to go, and only if the person asks
- `ytk-codon-optimise` to rewrite the gene for yeast, and only if the person asks
- `ytk-design-primers` to amplify it from cDNA
- `ytk-add-overhangs` to order the flanked fragment instead
