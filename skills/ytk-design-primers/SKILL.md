---
name: ytk-design-primers
description: Design PCR primers that add the MoClo Yeast Toolkit (YTK) Golden Gate flanks to a coding sequence, so the product amplified from cDNA is already YTK-compatible and drops straight into the pYTK001 entry vector. Works out the binding region, Tm, GC and the annealing temperature, and writes one CSV to order from. Use this whenever someone wants primers, oligos, a primer pair, something to order from a supplier, or asks how to amplify a gene for YTK cloning, Golden Gate or the entry vector - including when they only say "design primers for these genes", "what oligos do I order", or "make me primers for this batch". For the whole job at once, use ytk-workflow.
---

# Design YTK primers

Takes one coding sequence or a batch of them. Gives back one CSV of primer
pairs, ready to order.

The primer carries the YTK flanks on its 5' end, so the PCR product is already
YTK-compatible. Amplification is from cDNA.

## What you need from the person

- **The sequences.** A file, or one sequence on the command line.
- Nothing else. This skill only makes Type 3 parts, so there is no part type to
  choose.

Never guess a sequence, and never change one.

## Run it

```
python3 "$CLAUDE_SKILL_DIR/scripts/design_primers.py" --input genes.csv
```

One sequence instead of a file: write it to a one-line CSV first, so there is
only ever one input path to reason about.

```
printf 'name,sequence\nmy_gene,ATGAAA...TAA\n' > genes.csv
```

Input can be `.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`.

**Map files:** a plasmid map holds the insert plus the vector's own machinery,
and a resistance marker is usually labelled `CDS` while the insert is labelled
`misc_feature`. So markers and origins are set aside before the insert is
picked: first by matching a feature's sequence against the published parts in
`data/ytk_parts.tsv`, then by name against `data/backbone_features.tsv`. If one
feature is left, it is used. If several are, the run stops and lists them, and
`--feature NAME` says which to take.

## What comes out

`out/ytk_primers.csv`, two rows per sequence, one per oligo, named
`<name>_forward` and `<name>_reverse`. The columns:

`oligo`, `sequence`, `full_length`, `binding_region_length`,
`bind_region_gc`, `tm_phusion`, `warnings`.

**Every number except `full_length` is measured on the binding region alone**,
not the whole oligo. Only the binding region sticks to the template in the first
PCR cycle, and the flanks are fixed YTK adapters, so their GC is not the
designer's to fix.

**There is no combined Tm column.** An annealing temperature belongs to a pair,
and a row is one oligo, so the column was only ever right when the oligos
happened to be ordered two by two. The suggested temperature - three degrees
above the lower of the two primer Tms, which is NEB's rule for primers over 20
nt - is printed on screen instead, per pair. Confirm it on the NEB calculator
before you order.

`warnings` holds **that oligo's** problems only - a low Tm, a GC outside the
band, a long single-base run, an over-length primer. A run of identical bases in
the forward primer says nothing about the reverse one, so it is not repeated
there. The one exception is an unbalanced pair, which is written to both rows
because it is a fact about both. A cut site inside the gene is **not** here:
that is a fact about the gene, and it belongs to `ytk-cds-qc` and to
`fragments/summary.csv`.

## The rules it follows

| Rule | Value |
|---|---|
| Binding region | 20-25 bp aimed for, never below 18 |
| Total primer length | 50 bp target |
| Tm | 60-68 C aimed for, hard floor 50 C |
| GC | 40-60% ideal, 30-40% is fine |
| 3' end | must be C or G |
| Forward and reverse Tm | within 2 C of each other |

Runs of four or more identical bases are reported but never refused: a quarter
of the lab's working primers have one.

## The pad

The few spare bases at the outer end exist only so BsmBI has room to cut near
the end of a linear fragment.

- A primer uses **4**, and may go up to 10. Never fewer than 4.
- `ytk-add-overhangs` keeps the full **10** for a fragment ordered from a
  synthesis company.

**Those two outputs differ on purpose, and it does not change the plasmid.** The
pad sits outside both BsmBI sites, so it is cut off and thrown away. A 10-base
pad and a 4-base pad give the same part plasmid.

## Stop and ask

Unless the person already said what to do, stop and ask when:

- the sequence is not a valid Type 3 CDS - no ATG, no stop codon, or a length
  that is not a whole number of codons
- the sequence holds a BsmBI or BsaI site
- a primer would have to be longer than 50 bp
- no binding region at all can be found

In a batch, do not stop at the first problem. The script collects every problem,
prints one table, and asks once. Show the person that table.

**A pair that will not balance is not a refusal.** If nothing gets the two Tms
within 2 C, the closest pair still comes out, with a warning saying so. Use it as a
starting point and finish the design by hand. Point out the warning; do not present those primers as finished.

## Where the Tm comes from

primer3, with Phusion's reaction conditions, plus a fixed offset that was
measured once against the 45 sound rows of a 49-primer sheet, checked on the
NEB Phusion calculator.
Do not tune it per primer. `tests/test_primers.py` fails if it drifts.

## After this

The PCR product goes into pYTK001 with `ytk-clone`. If the PCR fails in the lab
after a few tries, order the flanked sequence instead with
`ytk-add-overhangs`.
