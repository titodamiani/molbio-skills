---
name: ytk-cds-qc
description: Check a coding sequence, or a batch of them, before any MoClo Yeast Toolkit (YTK) cloning is designed. Confirms each one is a valid Type 3 CDS - starts with ATG, ends with a stop codon, a whole number of codons - and finds any BsmBI or BsaI site inside that would break the Golden Gate reaction. Use this whenever someone wants a sequence checked, screened or validated before cloning, asks whether a gene is suitable for YTK or Golden Gate, asks if a CDS has internal cut sites or restriction sites, or says something like "are these genes ok", "check these sequences first", or "will this clone".
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
python3 skills/ytk-cds-qc/scripts/cds_qc.py --input genes.csv
```

One sequence instead of a file:

```
python3 skills/ytk-cds-qc/scripts/cds_qc.py --sequence ATGAAA...TAA --name my_gene
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

**Nothing, unless the person asks.** Report it and stop.

Removing a site means changing the DNA. Never do that on your own. If they ask
for it, the change must be silent - a different codon for the same amino acid -
and every changed base has to be reported.

## A stop codon in the middle

Reported as a warning, not a failure. It usually means one of two things: the
reading frame is wrong, or the sequence is not the coding sequence it was
thought to be. Worth saying out loud before anyone orders anything.

## After this

- `ytk-design-primers` to amplify it from cDNA
- `ytk-add-overhangs` to order the flanked fragment instead
