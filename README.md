# molbio-skills

Makes SnapGene plasmid maps for MoClo Yeast Toolkit cloning, so you do not
have to build them by hand.

You give it a list of genes. It puts the Golden Gate overhangs on each one,
works out what each finished plasmid looks like, writes a SnapGene `.dna` file
for every fragment and every plasmid, puts the part labels on each map, and
checks the results.

It only knows the Yeast Toolkit. That is on purpose.

## Install

In Claude Code, run these two lines:

```
/plugin marketplace add titodamiani/molbio-skills
/plugin install molbio-skills
```

That is the whole setup. You do not need to install Python packages. The
first time a skill runs it checks for what it needs and installs it.

## Use it

Just ask. For example:

> Add Type 3 overhangs to the genes in `candidates.fasta` and put the
> fragments in `out/`.

Or, from genes to finished maps:

> Take the genes in `candidates.fasta` and make me plasmid maps in `out/`.

Or, for the whole job at once:

> Run the full YTK job on `genes.csv` and check the results against my files
> in `refs/`.

Always say the part type — type 3 for a coding sequence, for example. The
tool never guesses it, because a wrong type puts the part in the wrong slot of
an assembly and nothing downstream would notice.

## What you give it

One file of genes. A CSV is best:

```
name,sequence,plasmid
Pi_fim_NCS_c1,ATGATTCCT...,pTP412
Pi_fim_NCS_c3,ATGGTTGCC...,pTP768
Pi_fim_OMT_c1,ATGGTCTTA...,pTP002
```

- **name** — what the gene is called
- **sequence** — the DNA
- **plasmid** — what to call the finished plasmid. Optional.

The plasmid column exists because real plasmid numbers jump about. pTP412,
pTP768 and pTP002 are fine side by side. The tool never makes up a number.

Leave the plasmid column out, or leave a cell blank, and that plasmid is
called `<gene>_<backbone>.dna` instead — for example
`Pi_fim_NCS_c1_pYTK001.dna`. The vector goes in the name so the file still
makes sense months later.

A header row is optional. Capitals in the header do not matter.

**FASTA** works too — `.fa` or `.fasta` — but a FASTA file has nowhere to put
a plasmid name, so you always get `<gene>_<backbone>.dna`.

Nothing else is accepted. If you hand it something else it says so rather
than guessing.

If two rows share a gene name, or share a plasmid name, it stops. Otherwise
one file would quietly overwrite the other.

Your sequences are never changed. If one looks wrong, the tool tells you and
stops.

## What you get back

From the overhangs step, one file per gene:

- `<gene>_oh.dna` — the fragment to order, with the flanks on it

Then, for each fragment, two more:

- `<fragment>.dna` — the fragment on its own, linear
- `<plasmid>.dna` — the finished circular plasmid, with labels

Open them in SnapGene as usual.

## The five skills

| Skill | What it does |
|---|---|
| `ytk-add-overhangs` | Puts the YTK flanks on a gene, ready to order. |
| `ytk-clone` | Builds the plasmid maps from flanked fragments. |
| `ytk-annotate` | Puts the part labels on a map. |
| `ytk-qc` | Checks a map, on its own or against a reference. |
| `ytk-batch` | Runs all four in order. |

Each one works on its own. You can also ask for just one of them.

`ytk-clone` takes the fragment as ordered, with the flanks already on it. Hand
it a bare gene and it stops and sends you to `ytk-add-overhangs`. It will not
add the flanks for you, because that would assume a fragment design you may
never have ordered.

## Two things worth knowing

**Where the map starts.** A plasmid is a loop, so it has no natural first
base. Every map this tool makes starts at the same point in the backbone.
That means the gene is never cut in half by the start of the file, and the
same gene always gives you the same map.

**Comparing maps.** Because a plasmid is a loop, two files can hold exactly
the same DNA while starting at different points. Compared as plain text they
look different, but they are not. `ytk-qc` turns one until it lines up, and
tells you `matches the reference, turned by N bp`. That is a pass, not a
failure.

**Genes with a cut site inside.** Some genes hold a BsmBI or BsaI site in the
coding sequence. The tool handles it and leaves the gene alone. It flags
those genes, because at the bench they cannot be re-cut with the same enzyme.

**Stop codons.** A coding sequence has to bring its own stop codon. The flanks
do not add one. If a gene has no stop codon the run stops, and you say
`--no-stop-codon <name>` if that was on purpose.

## The two tables

`data/ytk_overhangs.tsv` holds the flanks for every part type, one row each.
`ytk-add-overhangs` reads it. A row can also carry a note in the `unsure`
column, and then that part type stops the run until someone approves it. When
the overhangs are confirmed, clearing the cell turns the gate off.

`data/ytk_parts.tsv` holds every toolkit part: its name, its type, its
sequence and its two junctions. Both `ytk-clone` and `ytk-annotate` read it.

It was built from the published toolkit files with
`data/build_parts_table.py`. To rebuild it:

```bash
pip install pandas xlrd
python3 data/build_parts_table.py ~/Downloads/ytk/plasmids
```

Labels always come from this table, never from another map file. Map files in
circulation carry mistakes — one `pYTK001.dna` says the ColE1 origin covers
the whole plasmid, when it is 764 bp.

## Source

Lee ME, DeLoache WC, Cervantes B, Dueber JE.
*A Highly Characterized Yeast Toolkit for Modular, Multipart Assembly.*
ACS Synthetic Biology 2015, 4(9), 975–986.
[doi:10.1021/sb500366v](https://doi.org/10.1021/sb500366v)

Plasmids: Addgene kit #1000000061.

## Tests

```bash
python3 tests/test_ytk.py
```

The reference `.dna` files in `tests/data/` were made by hand in SnapGene.
The five cases each cover something different: a map that starts elsewhere on
the circle, a gene with a BsmBI site inside it, a gene with a BsaI site
inside it, a long insert, and the linear writer.

The genes in `tests/data/` are bare, and the tests flank them with
`ytk-add-overhangs` before cloning. So the two skills are tested together, and
the finished maps still have to match the files made by hand.

## Licence

MIT.
