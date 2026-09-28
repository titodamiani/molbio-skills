# ytk-skills

Does the MoClo Yeast Toolkit cloning work you would otherwise do by hand.

It works on one coding sequence or on a batch of them. It checks that each one is
a valid Type 3 coding sequence and has no internal BsmBI or BsaI site, designs the
PCR primers, simulates cloning the flanked sequence into a backbone, and writes
labelled plasmid maps, as SnapGene `.dna` or GenBank.

## Install

In Claude Code:

```
/plugin marketplace add titodamiani/ytk-skills
/plugin install ytk-skills
```

Then install the Python packages once:

```
python3 -m pip install -r requirements.txt
```

The skills do not install anything themselves. If something is missing they print
that exact line and stop, so the pinned versions are the versions you get.

## Scope: part Type 3 only

One part type is supported: **Type 3**, a whole coding sequence that starts with
`ATG` and ends with a stop codon.

**The part type is never guessed from a sequence.** Type 3 versus 3a versus 3b is
your design decision, not a property of the DNA. Nothing here will pick for you,
because a wrong type puts the part in the wrong slot of an assembly and nothing
downstream would notice.

## No tags on these parts

Read this before you build a set of these parts.

A part built here keeps the CDS's own stop codon, so:

| You want | With these parts |
|---|---|
| the plain protein | **works** |
| a C-terminal tag | **not possible** - translation stops before the tag |
| an N-terminal tag | **not possible** - the part fills the whole Type 3 slot |

The published Type 3 parts drop the stop codon and add `GG`, so a tag can be
fused on later and the stop comes from the next part. This plugin does not do
that, because it would mean editing your sequence.

Nothing you have already made is wrong. For promoter + CDS + terminator these
parts are correct. If you ever want tags, `NOTES.md` has the recipe, and it needs
no sequence editing either - you just supply a CDS without a stop codon.

## The seven skills

Each one takes one sequence or a batch.

**`ytk-cds-qc`** - are these valid Type 3 coding sequences, and do any hide a
BsmBI or BsaI site?

> Check the genes in `candidates.fasta` before I order anything.

> That gene has a BsmBI site. Remove it with a silent codon change and show me
> what you changed.

**`ytk-add-overhangs`** - puts the YTK flanks on a sequence, so you can order the
whole fragment from a synthesis company.

> Add Type 3 overhangs to the genes in `candidates.fasta` and put the fragments
> in `out/`.

**`ytk-design-primers`** - PCR primers that carry the flanks, so the product
amplified from cDNA is already YTK-compatible.

> Design primers for the genes in `PiperGenes.csv`. I'm amplifying from cDNA.

**`ytk-clone`** - cuts and joins a flanked fragment into a backbone and writes the
map. This is what SnapGene's cloning simulation does.

> Clone the fragments in `out/ytk_order.csv` into pYTK001 and give me the maps.

> Clone these into `pRS416.gb` with BsaI instead.

**`ytk-annotate-map`** - puts the part labels on a map file, `.dna` or GenBank.

> This plasmid map opens blank in SnapGene. Can you label it?

**`ytk-verify-map`** - checks a map file, on its own or against a reference.

> Check the maps in `out/maps/` against my hand-made ones in `refs/`.

**`ytk-workflow`** - all six, in order.

> Run the full YTK job on `genes.csv` and put everything in `out/`.

## The real lab workflow

1. Pick a coding sequence, or a batch of them.
2. **Check for BsmBI and BsaI sites inside.** Those break Golden Gate.
3. Design primers that carry the YTK flanks.
4. PCR from cDNA.
5. **If the PCR keeps failing**, order the flanked sequence from a synthesis
   company instead. That is what `ytk-add-overhangs` gives you.
6. Golden Gate into pYTK001 with BsmBI, and get a map.

## The pad, and why two outputs differ

At the outer end of every flank sit a few spare bases. They exist only so BsmBI
has room to cut near the end of a linear fragment.

- **`ytk-add-overhangs` keeps all 10.** That is the fragment you order.
- **`ytk-design-primers` trims to 4.** Shorter primers cost less.

**This is not a bug. Both give the identical plasmid.** The pad sits *outside*
both BsmBI sites, so it is cut off and thrown away. A fragment with a 10-base pad
and a PCR product with a 4-base pad give the same part. There is a test for it.

The pad is never below 4 bases. BsmBI needs that much room to cut.

## What you give it

Any of `.fa`, `.fasta`, `.csv`, `.gb`, `.gbk` or `.dna`. A CSV is easiest:

```
name,sequence,plasmid
Pi_fim_NCS_c1,ATGATTCCT...,pTP412
Pi_fim_NCS_c3,ATGGTTGCC...,pTP414
Pi_fim_OMT_c1,ATGGTCTTA...,pTP0457
```

- **name** - what the gene is called
- **sequence** - the DNA
- **plasmid** - what to call the finished plasmid. Optional, and carried all the
  way through, because real plasmid numbers are not in order.

For a `.gb` or `.dna` map, the insert is picked out for you: resistance markers
and origins are set aside first, by matching against the published parts and a
catalogue of common backbone features. If more than one candidate is left, it
stops and lists them, and you say which.

## What you get back

- `ytk_primers.csv` - two rows per gene, columns matching an oligo stock sheet
- `ytk_order.csv` - the flanked sequences to order
- `maps/` - one labelled circular map per plasmid
- `fragments/` - one linear map per fragment

Maps are SnapGene `.dna` by default. Add `--format genbank` to `ytk-add-overhangs`
and `ytk-clone` to get `.gb` files instead, for anyone who does not have SnapGene.
The plasmid is identical either way.

**Your sequences are never changed unless you ask.** There is one exception, and
you have to name it: `ytk-cds-qc --remove-sites` swaps a codon to remove an
internal BsmBI or BsaI site. It writes a new file, leaves yours alone, prints
every changed base, and refuses if the protein would change at all. On the two
real genes in `tests/data/` it changes one base.

## When it stops and asks

Not a valid Type 3 CDS. A BsmBI or BsaI site inside the sequence. A primer that
would have to be longer than 50 bp. A cloning whose sticky ends do not match.

**In a batch it never stops at the first problem.** Every skill collects the whole
batch's problems, prints one table, and asks once.

## The data tables

| File | Holds |
|---|---|
| `data/ytk_parts.tsv` | all 96 published parts, with sequences and junctions |
| `data/ytk_overhangs.tsv` | the adapter pair for each part type |
| `data/ytk_part_types.tsv` | the part-type overhangs, from the paper |
| `data/backbone_features.tsv` | markers and origins, so they are not mistaken for an insert |

Every overhang was checked twice: against the paper, and against all 96 published
plasmid maps in `reference/ytk_plasmids/` by finding their BsaI sites. Both agree
with every row. `data/ytk_parts.tsv` is generated from those maps, and a test
rebuilds it and compares byte for byte, so the table and the generator cannot
drift apart.

## Source

Lee MW, DeLoache WC, Cervantes B, Dueber JE (2015). *A Highly Characterized Yeast
Toolkit for Modular, Multipart Assembly.* ACS Synthetic Biology 4(9), 975-986.
[doi:10.1021/sb500366v](https://doi.org/10.1021/sb500366v)

Plasmids: Addgene kit #1000000061, all 96 committed in `reference/ytk_plasmids/`.
The two paper PDFs are not committed - 9 MB, and the 1 KB the code needs is in
`data/ytk_part_types.tsv`.

## Tests

```
python3 tests/test_ytk.py       # cloning, writing .dna and GenBank, annotation
python3 tests/test_primers.py   # the primer designer, vs 20 real primer pairs
python3 tests/test_parts.py     # the data tables and the entry vector
python3 tests/test_inputs.py    # the four formats, and asking once
```

`tests/test_primers.py` measures against a real oligo sheet and gene list. Those
hold unpublished sequences, so they are not in the repo: those tests are skipped,
with a message, when the files are not on your machine.

`NOTES.md` records the facts that are hard to find again.

## Licence

MIT.
