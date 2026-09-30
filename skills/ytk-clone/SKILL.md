---
name: ytk-clone
description: Simulate Golden Gate cloning and write the map of the finished construct, the way SnapGene's cloning simulation does. Cuts flanked fragments and a backbone with a Type IIS enzyme, checks the sticky ends match, joins them, and writes one map per plasmid, as a GenBank or SnapGene .dna file. Defaults to the MoClo Yeast Toolkit entry reaction - pYTK001 with BsmBI - and takes any other backbone and any Type IIS enzyme. Use this whenever someone has flanked fragments and wants plasmid maps, SnapGene files, Golden Gate assembly simulated, YTK Type 3 part plasmids, entry vector cloning, or asks whether a cloning will work - including when they only say "turn these into plasmids", "make maps for this batch", or "clone these into the vector". For bare genes with no overhangs yet, use ytk-add-overhangs first.
---

# Build YTK part plasmids

This skill turns flanked fragments into SnapGene plasmid maps.

The input is the fragment as ordered from a synthesis company, with the YTK
flanks already on it. A bare gene is refused. Use **ytk-add-overhangs** to
design the fragment first.

For each fragment it writes one map: `<plasmid>.gb`, the finished circular
plasmid.

## Run it

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/clone.py" --input FRAGMENTS --genes GENES
```

`FRAGMENTS` is `fragments/summary.csv` from **ytk-add-overhangs**. That clones
into pYTK001 with BsmBI, which is the YTK entry reaction.

`GENES` is the input file: `input.csv`, or
`corrected_genes/input_corrected.csv` when sites were removed. The plasmid names
and the record of what was done to each gene come from there, because both are
facts about the gene rather than about the fragment.

Without `--outdir` the maps go in a `ytk_output/` folder beside the input, and
the script prints where. An input that is already inside one goes back into it,
so chaining the steps by hand does not bury a folder in the last one.

Any other vector, any other Type IIS enzyme:

```bash
python3 "$CLAUDE_SKILL_DIR/scripts/clone.py" --input FRAGMENTS \
    --backbone my_vector.gb --enzyme BsaI
```

`--backbone` takes `.dna`, `.gb`, `.gbk` or FASTA. `--enzyme` takes any name
Biopython knows. Away from the default pair, the Type 3 junction check is
skipped, because `TATG`/`ATCC` only means something for the entry reaction. The
overhangs are still checked: pydna refuses to join ends that do not match.

`$CLAUDE_SKILL_DIR` is this skill's own folder, so the command works from any
directory. The script finds the shared code and the parts table by itself.

## What comes out

One folder, named after the backbone actually cloned into, so a run against
another vector cannot be mistaken for an entry-vector run:

```
plasmids/
  pYTK001/
    summary.csv    plasmid, sequence, part_type, codon_opt, notes
    <plasmid>.gb   one circular map per plasmid, not yet labelled
```

`sequence` is the gene name: which gene is in that plasmid. `part_type` comes
from the fragment file, and `codon_opt` and `notes` from `--genes`. All are
optional: a bare FASTA of fragments still clones, it just leaves them blank.

**The format is not chosen here.** It is read off the fragment maps beside the
input, so the maps come out in whatever format the fragments were written in.
`--format` overrides that if you really need to. The plasmid is identical either
way; only the file format changes.

So the next step can just take `plasmids/pYTK001/*.gb` from the output folder.

The linear post-digest pieces are not written. Nothing downstream reads them,
and the fragment maps from `ytk-add-overhangs` are the ones worth keeping.

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

Plasmid names come from the plasmid column. Real lab numbers are not in order.
pTP412, pTP768 and pTP002 can be three plasmids in a row. So the script never
makes up a number.

Where a row gives no plasmid name, and for every FASTA file, the plasmid is
named `<gene>_<backbone>.gb` — for example `Pi_fim_NCS_c1_pYTK001.gb`. The
vector is in the name because a folder of `_plasmid.gb` files does not say
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
unexpected stop — say so and stop. Do not fix it, even when the fix is small
and invisible, such as a codon change to remove a cut site. A corrected gene
looks right, so it reaches the bench unnoticed. Only the owner of the gene
decides.

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
starts the map at base 1 of the backbone.

Two things follow from that, and both matter:

- The start is always in the backbone, so the gene is never split across the
  start of the file. A split gene looks broken in SnapGene even when the DNA
  is fine.
- The same gene always gives the same map, so two people get the same file.

One catch with a backbone read straight from a kit: its base 1 can sit inside the
piece that drops out, and then it is not in the finished plasmid at all. When
that happens the map starts at the first base of the piece that was kept. The
plasmid is the same either way, just numbered from a different point.

## After this

The maps have no labels yet. Use **ytk-annotate-map** to label the parts, then
**ytk-verify-map** to check the results. **ytk-workflow** does the whole job in
order.
