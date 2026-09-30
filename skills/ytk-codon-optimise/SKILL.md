---
name: ytk-codon-optimise
description: Codon optimise a coding sequence for Saccharomyces cerevisiae, rewriting each codon as the most-used codon for its amino acid while keeping the protein exactly as it was and leaving the gene its own start and stop codon. Reports the yeast CAI, GC and the longest repeat before and after, and never touches the input file. Use this whenever someone wants a gene codon optimised, codon-optimised for yeast, recoded, "made to express in yeast", or asks for a yeast version of a bacterial, plant, insect or mammalian gene - including when they only say "optimise these genes", "codon optimise for S. cerevisiae", or "rewrite the codons". Rewriting codons can create a BsmBI or BsaI site, so ytk-remove-cut-sites must run after this, every time. To only remove a cut site and change as little as possible, use ytk-remove-cut-sites on its own. For the whole job at once, use ytk-workflow.
---

# Codon optimise a coding sequence for yeast

Rewrites every codon as the most-used codon for its amino acid in *S. cerevisiae*.
The protein is translated before and after and compared, and if it differs at all
nothing is kept. The gene keeps its own `ATG` and its own stop codon.

**If the gene is going to a synthesis company anyway, their optimiser is better
than this one and free.** Twist, IDT and GenScript all optimise at order time,
against their own synthesis model, and they will not build a sequence their own
screen rejects. Use this when the sequence is needed before the order goes in,
when exactly which codons moved has to be on the record, or when the gene must
come out of the plugin's own pipeline.

Run `ytk-cds-qc` first and show its table. That is where `CAI` comes from, and a
gene already at a high CAI has nothing to gain here.

## Run it

```
python3 "$CLAUDE_SKILL_DIR/scripts/codon_optimise.py" --input genes.csv
```

Input can be `.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`. Another organism's
table with `--taxid`, which takes any taxonomy id the `python_codon_tables`
package carries - `4932` is yeast and the default, `316407` is *E. coli*, `9606`
is human.

## Run ytk-remove-cut-sites after this, always

Rewriting codons can put a `BsmBI` or `BsaI` site into a gene that had none, and
nothing here looks for one. Removing a site is `ytk-remove-cut-sites`' job and
only its job, so there is one implementation of it in the plugin rather than two.

This script exits 1 when a site is present afterwards and prints the exact
command to run, so the step cannot be skipped by accident. Run it even when the
script says no site was created, so the check is on the record.

## Output

An `optimised_genes/` folder. The input file is never touched.

| File | What it is |
|---|---|
| `input_optimised.csv` | the handoff to the next step, in the same shape as `input.csv` |
| `summary.csv` | for reading: the old sequence beside the new one |
| `changes.csv` | one row per changed codon: the codon number, the amino acid, and what it was |

`input_optimised.csv` has the same columns as `input.csv` and as
`corrected_genes/input_corrected.csv`, on purpose: whichever file the next step
is pointed at, it reads the same columns, so nothing downstream has to know which
fixes ran. The measurements on it are retaken from the new sequence, so the row
describes the gene it is on and not the one it came from.

`codon_opt_method` holds the table used, for example
`argmax/4932/python_codon_tables-0.1.18`. The package version is part of it
because the table lives in the package: a bump that moved one codon would
otherwise change a sequence with nothing in the output saying why.

## Two things get worse, and neither is fixed

Both are printed. Nothing in this plugin fixes either.

- **GC falls, every time.** Yeast's preferred codons are AT-rich, so a gene at
  50% GC comes out near 32%. Low GC makes PCR and synthesis harder. A warning
  prints under 35%. On the four real genes in `tests/data/` it lands at 30 to 33%,
  so expect the warning rather than being surprised by it. **Say this to the
  person before they order anything.**
- **Repeats grow.** The same peptide motif twice becomes the same DNA twice, and
  a synthesis company can refuse the order over it. A run of one amino acid is
  handled - the same codon is never used more than twice running, so a poly-Q
  tract comes out `CAACAACAG` and not `CAACAACAA` - but a repeated motif is not.

## Never optimise a sequence unless the person asked in that message

This changes hundreds of bases. Never run it on your own initiative, never to get
a failing check to pass, and never because a gene "looks bacterial".

Show the person the `CAI` column from `ytk-cds-qc` and the number of changed
codons before saying the job is done. If `CAI` was already around 0.5 or higher,
say plainly that the rewrite probably buys nothing.

**One deviation from the house rule, on purpose.** `ytk-remove-cut-sites` prints
every changed base, because it moves one. A rewrite moves hundreds, so the
console prints the full table only for a gene with 40 or fewer changed codons,
and `changes.csv` always has every row. Point at that file rather than claiming
every base was shown.

## When not to optimise at all

- **The `CAI` is already high.** Nothing to gain.
- **The gene is from yeast, or from a close relative.** A *Kluyveromyces* or
  *Candida* gene scores respectably and would be rewritten for nothing. `CAI`
  cannot tell "already yeast" from "yeast-like".
- **Primers or oligos already exist for this gene.** They will not match the new
  sequence. Someone with oligos in the freezer now has oligos for a gene that no
  longer exists.
- **Anything overlapping the CDS matters.** An internal promoter, a uORF, an
  intron, a restriction site someone was relying on, a programmed frameshift: a
  rewrite destroys all of it silently, because it only protects the protein.
- **The reading frame is not certain.** A sequence in the wrong frame can still
  be a whole number of codons that starts `ATG` and ends in a stop, and it
  rewrites into something that looks clean and means nothing. A stop codon in the
  middle stops this script, but that only catches the wrong frames that happen to
  contain one.

## What it does not do

- **Not codon harmonisation.** That is a different algorithm - matching the
  source organism's rarity profile rather than the target's most-used codon - and
  this does not do it.
- No GC targeting, no repeat avoidance, no secondary-structure work. Both of the
  first two need a search over the whole sequence and a target nobody has set.
- No cut-site avoidance. That is the next step's job.
- No protein sequence as input. There has to be an input gene to compare against.
- Type 3 only. Rewriting codons needs a reading frame.

## After this

- `ytk-remove-cut-sites` - **not optional**
- then `ytk-add-overhangs`, reading `input_optimised.csv` and never the original
- `ytk-design-primers`, also reading the new file: primers designed from the old
  sequence would not bind the rewritten gene
