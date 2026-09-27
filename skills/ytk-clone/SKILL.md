---
name: ytk-clone
description: Build SnapGene .dna plasmid maps from flanked YTK fragments using the MoClo Yeast Toolkit (YTK). Reads a FASTA or CSV file of fragments that already carry the YTK overhangs, cuts and joins each one into the pYTK001 entry vector, and writes one .dna file per fragment and one per finished plasmid. Use this whenever someone has ordered or designed YTK fragments and wants plasmid maps, SnapGene files, Golden Gate assembly, YTK Type 3 part plasmids, or entry vector cloning - including when they only say something like "turn these into plasmids", "make maps for this batch", or "clone these into the toolkit vector". For bare genes with no overhangs yet, use ytk-add-overhangs first.
---

# Build YTK part plasmids

This skill turns flanked fragments into SnapGene plasmid maps.

The input is the fragment as ordered from a synthesis company, with the YTK
flanks already on it. A bare gene is refused. Use **ytk-add-overhangs** to
design the fragment first.

For each fragment it writes two files:

- `<fragment>.dna` — the fragment on its own, linear
- `<plasmid>.dna` — the finished circular plasmid

## Run it

```bash
python3 scripts/clone.py --input FRAGMENTS --outdir OUTDIR
```

`scripts/clone.py` sits in this skill's own folder. The script finds the
shared code and the parts table by itself, so it works wherever the plugin is
installed.

## Input

Only two kinds of file are accepted:

- **CSV**, ending in `.csv`, with a name column, a sequence column, and an
  optional plasmid name column
- **FASTA**, ending in `.fa` or `.fasta`

A CSV may have a header row or not. Header names are matched whatever the
capitals, so `Name`, `NAME` and `name` all work. The plasmid column may be
called `plasmid`, `plasmid_name` or `construct`.

```
name,sequence,plasmid
Pi_fim_NCS_c1,actcgacaacCGTCTCatcGGTCTCaTATGATTCCT...gttgtggtgt,pTP412
Pi_fim_NCS_c3,actcgacaacCGTCTCatcGGTCTCaTATGGTTGCC...gttgtggtgt,pTP768
Pi_fim_OMT_c1,actcgacaacCGTCTCatcGGTCTCaTATGGTCTTA...gttgtggtgt,pTP002
```

The sequence column holds the flanked fragment, not the bare gene. The names
stay the gene names, because the map is labelled and filed under those.

Anything else is refused with a message saying so. If you are handed a
GenBank file, an Excel sheet or a Word document, do not convert it quietly.
Tell the person what you have and ask what they want, because guessing at
someone's sequences is how errors get built into a plasmid.

## Plasmid names

Plasmid names come from the plasmid column. Real lab numbers are not a tidy
series — pTP412, pTP768 and pTP002 sit happily side by side — so the script
never counts them out for itself.

Where a row gives no plasmid name, and for every FASTA file, the plasmid is
named `<gene>_<backbone>.dna` — for example `Pi_fim_NCS_c1_pYTK001.dna`. The
vector is in the name because a folder of `_plasmid.dna` files does not say
what anything is in.

Prefer the CSV with a plasmid column. If someone hands you a FASTA and cares
what the plasmids are called, ask for the names rather than inventing them.
A made-up plasmid number that reaches a lab notebook is hard to undo.

Two rows sharing a gene name, or sharing a plasmid name, stop the run. One
file would otherwise be written over the other and nobody would see it
happen.

## Never change the sequences

Treat every input sequence as read-only.

The script enforces this itself, in three places:

1. It stops if a sequence holds anything other than A, C, G and T.
2. It checks the piece it cut out is in the fragment exactly as given, and
   then that it is in the finished plasmid. If either fails, the run stops.
3. It builds every plasmid before writing anything, so a failure leaves no
   files at all rather than half a batch.

Check 2 is the important one. A changed base leaves a plasmid that still
closes into a loop and still looks perfectly normal in SnapGene. Nothing
downstream would catch it.

Do not work around any of these. If you spot a problem with a sequence — a
missing start codon, a length that is not a whole number of codons, an
unexpected stop — say so and stop. Do not fix it, even when the fix is
obvious and silent, such as a codon change to remove a cut site. A corrected
gene looks right and ends up at the bench. Only the person whose gene it is
gets to decide.

## What the script does

1. Cuts the fragment with BsmBI and keeps the piece between the two designed
   cuts.
2. Checks that piece is a Type 3 part. See below.
3. Cuts the pYTK001 entry vector with BsmBI and keeps the larger piece.
4. Joins the two and closes the loop.
5. Turns the loop so it starts at base 1 of the backbone.

## What it refuses

**No pair of BsmBI sites.** The fragment looks like a bare gene. The script
says so and stops. It does not add the flanks for you, even if the person
says to go ahead — that would assume a fragment design they may never have
ordered. Send them to **ytk-add-overhangs** instead.

**The wrong part type.** The script reads the part's own overhangs and
compares them with the Type 3 pair, `TATG` and `ATCC`, taken from
`data/ytk_parts.tsv`. Anything else does not belong in this entry vector, so
the run stops and reports the overhangs it found.

Those overhangs come from the inner BsaI sites, not the outer BsmBI ones. The
BsmBI cut gives the same ends on every part type, because those ends are what
fits the entry vector.

## Fragments with an internal cut site

Some genes hold a BsmBI or BsaI site inside the coding sequence. Cutting then
gives extra pieces in the middle. The script joins those pieces back together,
so only the two designed outer cuts count and the gene comes through unchanged.

The gene is never edited to remove such a site. Cut sites are left exactly
where they are.

The output table has an `internal site` column, and the script names the
flagged genes at the end of the run. A `yes` is normal and needs no action.
Pass it on to the person, because at the bench that gene cannot be re-cut with
the same enzyme later.

## Where the map starts

A circle has no natural first base, so the file has to pick one. The script
always starts the map at base 1 of the backbone.

Two things follow from that, and both matter:

- The start is always in the backbone, so the gene is never split across the
  start of the file. A split gene looks broken in SnapGene even when the DNA
  is fine.
- The same gene always gives the same map, so two people get the same file.

## After this

The maps have no labels yet. Use the **ytk-annotate** skill to label the
parts, then the **ytk-qc** skill to check the results. The **ytk-batch** skill
does all three in order.
