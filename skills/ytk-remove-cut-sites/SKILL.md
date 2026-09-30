---
name: ytk-remove-cut-sites
description: Remove an internal BsmBI or BsaI site from a coding sequence by swapping one codon for a different codon of the same amino acid, so the protein is unchanged and the gene survives a MoClo Yeast Toolkit (YTK) Golden Gate reaction. Usually moves a single base. Every changed base is printed, the protein is compared before and after, and the input file is never touched. Use this whenever someone asks to remove, fix, get rid of or silently mutate an internal cut site, restriction site, BsmBI site or BsaI site in a gene, says something like "take the site out" or "make it Golden Gate compatible", or agrees when asked whether to remove one. It must also run straight after ytk-codon-optimise, every time, because rewriting codons can create a site. To find out whether a sequence has one at all, use ytk-cds-qc. For the whole job at once, use ytk-workflow.
---

# Remove a cut site from a coding sequence

A `BsmBI` or `BsaI` site inside a coding sequence gets cut during Golden Gate, so
the part comes apart. It is removed by swapping one codon for a different codon
of the same amino acid, which leaves the protein exactly as it was.

On the two real genes in `tests/data/` this needs **one base each**, both at the
third position of a codon where the genetic code is redundant.

## Run it

```
python3 "$CLAUDE_SKILL_DIR/scripts/remove_cut_sites.py" --input genes.csv
```

Input can be `.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`. For a map file,
`--feature NAME` says which feature holds the CDS.

## Run this after ytk-codon-optimise, always

Rewriting a gene's codons can create a site that was not there before, and
`ytk-codon-optimise` does not look for one. It exits 1 when a site is present and
prints the command to run. Run this step even when it says no site was created,
so the check is on the record.

Anything the earlier step did travels on the row with the sequence, so this
carries `codon_opt_method` forward rather than dropping it. Point it at
`optimised_genes/input_optimised.csv`, never at the original input.

## Safeguards. None of them is optional

- **The input file is never touched.** A `corrected_genes/` folder is written in
  the output folder, with two files:
  - `summary.csv` - `name`, `input_sequence`, `new_sequence`, `codon_opt`,
    `notes`. The old sequence beside the new one, for reading.
  - `input_corrected.csv` - the handoff to the next step, in the same shape as
    `input.csv`, so nothing downstream has to know this step ran. CSV and not
    FASTA, because a FASTA header cannot carry any of the extra columns.
- **The facts travel in their own columns.** `removed` holds the enzyme names, so
  a later step rebuilds the note rather than reading the sentence apart. Which
  enzymes were cleared cannot be measured after the swap, so that one has to be
  carried; what is still in the gene is measured again at every step.
- **The note says what was done and what is left** - for example
  `BsmBI site removed`, or `BsmBI site removed; BsaI site in the CDS` when no
  synonymous codon removes the second one. It is empty only when the gene was
  clean and nothing was done.
- **Every changed base is printed** with its position, its codon number, the old
  and new base, and the amino acid.
- **The protein is translated before and after and compared.** If the protein
  changes at all, nothing is kept and the run stops.
- **The start codon and the stop codon are never swapped.**
- If no synonymous codon removes a site, it says so, leaves that sequence alone,
  and exits 1 with the gene named. Those need doing by hand.

Show the person the change table.

## Never change a sequence to make a check pass

This skill exists because removing a cut site is sometimes the right call in the
lab. It is never the right call to silence a warning. **Never run it on your own
initiative** - only when the person asks in that message, or when
`ytk-codon-optimise` has just run.

`ytk-cds-qc` reports a site and cannot change anything. That is the right place
to stop and ask.

## Not a valid Type 3 CDS

Nothing is changed, and the run stops with the list. Swapping a codon needs a
reading frame, so a sequence has to be a valid coding sequence before this is
safe. Fix the sequence, or find out it is not the CDS anyone thought it was.

## After this

- `ytk-add-overhangs`, reading `input_corrected.csv`
- `ytk-design-primers`, also reading `input_corrected.csv`: a primer designed
  from the old sequence would carry the base that was just changed
