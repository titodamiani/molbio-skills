---
name: ytk-design-primers
description: Design PCR primers that add the MoClo Yeast Toolkit (YTK) Golden Gate flanks to a coding sequence, so the product amplified from cDNA is already YTK-compatible and drops straight into the pYTK001 entry vector. Works out the binding region, Tm, GC and the annealing temperature, and writes one CSV to order from. Use this whenever someone wants primers, oligos, a primer pair, something to order from a supplier, or asks how to amplify a gene for YTK cloning, Golden Gate or the entry vector - including when they only say "design primers for these genes", "what oligos do I order", or "make me primers for this batch".
---

# Design YTK primers

Takes one coding sequence or a batch of them. Gives back one CSV of primer
pairs, ready to order.

The primer carries the YTK flanks on its 5' end, so the PCR product is already
YTK-compatible. Amplification is from cDNA.

## What you must be told

- **The sequences.** A file, or one sequence on the command line.
- Nothing else. This skill only makes Type 3 parts, so there is no part type to
  choose.

Never guess a sequence, and never change one.

## Run it

```
python3 skills/ytk-design-primers/scripts/design_primers.py \
    --input genes.csv --outdir out/
```

One sequence instead of a file:

```
python3 skills/ytk-design-primers/scripts/design_primers.py \
    --sequence ATGAAA...TAA --name my_gene --outdir out/
```

Input can be `.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`. A map file with
several features will stop and list them; pass `--feature NAME` to say which
one holds the coding sequence.

## What comes out

`out/ytk_primers.csv`, two rows per sequence, named `<name>_forward` and
`<name>_reverse`. The columns match the lab's oligo stock sheet:

`name`, `sequence`, `Full length (bp)`, `Binding region length (bp)`, `GC%`,
`Tm Phusion (C)`, `Combined Tm Phusion (C)`, `warnings`.

**Every number is measured on the binding region alone**, not the whole primer.
Only the binding region sticks to the template in the first PCR cycle.

`Combined Tm Phusion (C)` is not a melting temperature. It is the annealing
temperature to run the PCR at: three degrees above the lower of the two primer
Tms, which is NEB's rule for primers over 20 nt.

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
- no primer pair meets the rules
- a primer would have to be longer than 50 bp

In a batch, do not stop at the first problem. The script collects every problem,
prints one table, and asks once. Show the person that table.

## Where the Tm comes from

primer3, with Phusion's reaction conditions, plus a fixed offset that was
measured once against 49 real primers checked on the NEB Phusion calculator.
Do not tune it per primer. `tests/test_primers.py` fails if it drifts.

## After this

The PCR product goes into pYTK001 with `ytk-clone`. If the PCR fails in the lab
after a few tries, order the flanked sequence instead with
`ytk-add-overhangs`.
