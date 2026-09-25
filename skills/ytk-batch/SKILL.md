---
name: ytk-batch
description: Run the whole MoClo Yeast Toolkit (YTK) job for a list of genes in one go - build the plasmid maps, label the parts, and check the results. Use this whenever someone has a FASTA or CSV file of genes and wants finished, labelled, checked SnapGene .dna files without running each step themselves, or says something like "do the whole thing", "the full run", "all of it", "start to finish", or "just give me the maps for these genes".
---

# Run the whole job

This skill runs the three other YTK skills in order for a batch of genes.

It has no scripts of its own. It uses the other skills, so each one keeps
working on its own if someone installs only that one.

## The order

**1. Build the maps** — use the **ytk-clone** skill.

Read its SKILL.md and run its script with the gene file and an output folder.
That writes one linear `.dna` per gene and one circular `.dna` per plasmid.

**2. Label the maps** — use the **ytk-annotate** skill.

Read its SKILL.md and run its script over the plasmid files from step 1, with
`--genes` pointing at the same gene file. That puts the part labels and the
gene name on each map.

**3. Check the results** — use the **ytk-qc** skill.

Read its SKILL.md and run its script over the files from step 2, with
`--genes` pointing at the same gene file. If the person has reference files
made by hand, run it once per file with `--reference` as well.

## Before you start

Ask for two things if they are not already clear:

- the gene file, CSV or FASTA
- where the output should go

A CSV can carry a plasmid name column, and that is the better input. If you
are given a FASTA, or a CSV with no plasmid column, the plasmids get named
`<gene>_<backbone>`, such as `Pi_fim_NCS_c1_pYTK001`. Say so, and ask whether
that is what they want.

Never invent plasmid numbers. Real lab numbers jump about, they are tracked
somewhere outside this tool, and a made-up one that reaches a notebook is
hard to undo. Ask instead.

## Stop on a problem

If a step fails, stop and report it. Do not carry on to the next step with a
broken file, and do not work around a failure by changing anyone's sequence.
A map that quietly went wrong costs far more later, at the bench, than a run
that stopped early.

## Report at the end

Say plainly:

- how many genes were read
- how many files were written, and where
- which genes hold an internal BsmBI or BsaI site (normal, but worth knowing,
  because those genes cannot be re-cut with the same enzyme later)
- the result of every check, and any that failed
