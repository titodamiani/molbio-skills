---
name: ytk-cds-qc
description: Check DNA sequences - one coding sequence or a batch of them - before any MoClo Yeast Toolkit (YTK) cloning is designed. Confirms each one is a valid Type 3 CDS - starts with ATG, ends with a stop codon, a whole number of codons - and finds any BsmBI or BsaI site inside that would break the Golden Gate reaction. Use this whenever someone wants a sequence checked, screened or validated before cloning, asks whether a gene is suitable for YTK or Golden Gate, asks if a CDS has internal cut sites or restriction sites, or says something like "are these genes ok", "check these sequences first", or "will this clone". For a finished plasmid map file use ytk-verify-map instead. For the whole job at once, use ytk-workflow.
---

# Check a coding sequence

Reads sequences and answers two questions:

1. **Is it a valid Type 3 CDS?** Starts with `ATG`, ends with a stop codon, and
   a whole number of codons.
2. **Does it hide a cut site?** A `BsmBI` or `BsaI` site inside the sequence
   gets cut during Golden Gate, so the part comes apart.

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

Nothing is written. Nothing is changed. Exits 1 if anything failed, so it can
gate a longer run.

## What it reports

A table, one row per sequence, then one table of everything wrong.

| Verdict | Meaning |
|---|---|
| `ok` | a clean Type 3 CDS with no internal cut site |
| `usable, with warnings` | a valid CDS, but something needs saying |
| `not a Type 3 CDS` | cannot be a Type 3 part as it stands |

All the faults are listed, not just the first, so a sequence can be fixed in one
go.

## What to do about a cut site

**Report it and stop. Do not remove it unless the person asks in that message.**

If they do ask, there is a flag for it:

```
python3 "$CLAUDE_SKILL_DIR/scripts/cds_qc.py" \
    --input genes.csv --remove-sites
```

It swaps one codon for a different codon of the same amino acid, so the protein
is unchanged. On the two real genes in `tests/data/` this needs **one base each**,
both at the third position of a codon where the genetic code is redundant.

Safeguards. None of them is optional:

- **The input file is never touched.** A new `corrected_genes/` folder is
  written in the output folder, with two files:
  - `summary.csv` - `name`, `input_sequence`, `new_sequence`, `codon_opt`,
    `notes`. The old sequence beside the new one, for reading.
  - `input_corrected.csv` - `name`, `plasmid`, `sequence`, `removed`,
    `codon_opt_method`. The handoff to the next step. CSV and not FASTA,
    because a FASTA header cannot carry any of the extra columns.
- **The facts travel in their own columns.** `removed` holds the enzyme names,
  so a later step rebuilds the note rather than reading the sentence apart.
  Which enzymes were cleared cannot be measured after the swap, so that one has
  to be carried; what is still in the gene is measured again at every step.
- **The note says what was done and what is left** - for example
  `BsmBI site removed`, or `BsmBI site removed; BsaI site in the CDS` when no
  synonymous codon removes the second one. It is empty only when the gene was
  clean and nothing was done.
- **Every changed base is printed** with its position, its codon number, the old
  and new base, and the amino acid.
- **The protein is translated before and after and compared.** If the protein
  changes at all, nothing is kept and the run stops.
- **The start codon and the stop codon are never swapped.**
- If no synonymous codon removes a site, it says so, and leaves the sequence
  alone.

Show the person the change table. Never apply this flag on your own initiative,
and never to get past a failing check.

## Never change a sequence to make a check pass

The flag above exists because removing a cut site is sometimes the right call in
the lab. It is never the right call to silence a warning. If a sequence fails,
report it.

## A stop codon in the middle

Reported as a warning, not a failure. It usually means one of two things: the
reading frame is wrong, or the sequence is not the coding sequence it was
thought to be. Say this before anyone orders DNA.

## After this

- `ytk-design-primers` to amplify it from cDNA
- `ytk-add-overhangs` to order the flanked fragment instead
