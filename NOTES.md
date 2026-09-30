# Notes for whoever works on this next

Things that are true, are not obvious, and would cost you a day to rediscover.
None of them is a bug in this repo.

## Type 3 is the default, not the limit

`ytk-add-overhangs`, `ytk-design-primers` and `ytk-clone` all take `--type`, and
types 1 to 8b work. Type 3 is the default: a whole coding sequence that starts
with `ATG` and ends with a stop codon.

The part type is never guessed from a sequence, because Type 3 versus 3a versus
3b is a design decision, not a property of the DNA. Every script prints the type
it used, so a wrong default cannot pass unnoticed, and `ytk-clone` checks the
fragment's own overhangs against the declared type on every run.

`data/ytk_overhangs.tsv` holds the adapters for every type, and
`data/ytk_part_types.tsv` holds the published overhangs. `lib/flanks.py` reads
both. Junctions come from `ytk_part_types.tsv` and never from the adapters: a
part type is defined by those four bases on each side and by nothing else, so
deriving them from a sequence would only ever agree with itself.

`data/ytk_parts.tsv` is not a junction source. It mixes in unpublished types
(`234`, `234r`, `678`, `cassette`), and the first row matching type 3 is a gene,
not a junction definition. That mistake once made the check compare a part
against `mTurquoise2`.

**`ytk-cds-qc` stays Type 3 only, and that is deliberate.** Its whole identity is
the question *is this a whole coding sequence*. `lib/silent.py` in particular
cannot be type-generic: finding a synonymous codon needs a reading frame. For any
other part type, skip that step and take the internal-cut-site information from
the `notes` column that `ytk-add-overhangs` writes.

**`custom` has `NNNN` adapters.** `flanks.junctions("custom")` stops the run, and
`lib/snapgene.py` refuses any sequence with a letter outside ACGT, so a custom
fragment could never be read back anyway. The row stays because it records the
paper's advice. Finish a custom part by hand.

## Getting a C-terminal tag later

Parts built here **cannot take a C-terminal tag**: the CDS keeps its own stop
codon, so translation halts before anything downstream. An N-terminal tag is out
too, because the part fills the whole Type 3 slot.

To add fusions, no code has to edit any DNA. Both recipes start the same way:
supply a CDS **without** a stop codon. They differ in which slot the tag takes.

**Route A: the tag is a Type 4a part.** The CDS stays a Type 3 part.

1. Use `ggATCC` as the right adapter, not `ATCC`. The `gg` is a frame filler, not
   a stop: `ATCC` is four bases and would break the reading frame, and `GGA TCC`
   is Gly-Ser, which is also a BamHI site. This is what the paper asks for, and
   what all five published Type 3 parts do.
2. The stop then comes from the next part. A **Type 4** terminator stops
   translation at once. A **Type 4a** tag (pYTK057-060: mTurquoise2, Venus,
   mRuby2, 6xHis-3xFlag) adds the tag and carries the stop at its own end, and
   then needs a **Type 4b** terminator after it.

**Route B: the tag is a Type 3b part.** The CDS becomes a **Type 3a** part
(`TATG` upstream, `ggTTCT` downstream). The tag sits in the second half of the
same slot and carries the stop codon. The published Type 3b tags are
**pYTK044-046** (mTurquoise2, Venus, mRuby2). A plain Type 4 terminator follows.

Route B has fewer tags to pick from. There are three Type 3b parts, and all
three are fluorescent proteins. There is no 6xHis-3xFlag among them: pYTK040
carries that tag as a Type **3a** part, for the N-terminal side. So a C-terminal
purification tag has to come from Route A's pYTK060.

Route B is the exact mirror of the N-terminal recipe, where the CDS becomes a
**Type 3b** part (`TTCT` upstream, `ATCC` downstream) and the tag is Type 3a.

Both routes cost the same PCR: amplify the CDS off an existing construct with a
reverse primer that drops the stop codon and carries the new right adapter.
Route B keeps the tag inside the Type 3 slot, so the terminator choice stays
free. Route A keeps the CDS as a plain Type 3 part, so it still drops into any
assembly that has no tag. `data/ytk_overhangs.tsv` already holds the adapters for
3a, 3b, 4a and 4b. No code reads them.

## Plasmid_Generator.xlsx has two faults

The spreadsheet in `reference/` is where the flank strings came from. Two things
in it are wrong, and the repo does not copy either.

**1. An extra stop codon.** Cells `J25` (Type 3) and `J27` (Type 3b) store the
right piece as `TAGATCC`. The `TAG` is a second stop codon, added whether or not
the CDS already has one. `data/ytk_overhangs.tsv` stores plain `ATCC`.

Checked against reality: all 24 matched reverse primers in `oligo_stocks` end at
their gene's true 3' end with its own stop codon, and none carries the extra TAG.

**Do not test for this by searching a fragment for `TAGATCC`.** A CDS whose own
stop codon is `TAG` produces that string legitimately. Test the adapter instead.

**2. The Left and Right formulas skip two columns.** They are `CONCAT(C:F)` and
`CONCAT(J:M)`, which miss columns `G` and `I`, where extra adapters are stored.
So the advanced path silently drops:

| Type | Cell | Lost | Consequence |
|---|---|---|---|
| 2 | `I24` | `agatc` | loses the BglII site |
| 3a | `I26` | `gg` | **shifts the reading frame of the fusion** |
| 4 | `G28` | `taactcgag` | loses the stop codon and XhoI site |
| 4a | `I29` | `taactcgag` | same |

Rows 34 and 37 (Types 8 and 8a) were fixed to `CONCAT(C:G)` and `CONCAT(I:M)`,
which is why their NotI sites survive. **Type 3 is unaffected**, because it uses
the simple path, `C3 = CONCAT(C25:M25)`, which takes the whole row.

`data/ytk_overhangs.tsv` keeps all four adapters, so the repo is more correct
than the sheet. The sheet itself has been left alone.

## The oligo stock sheet has six arithmetic slips

`oligo_stocks - Sheet1.csv` states a binding-region length that does not match
its own sequence in four rows, and a full length in two:

| Item | Column | Stated | Actual |
|---|---|---|---|
| TD014R | Binding region | 25 | 24 |
| TD020R | Binding region | 19 | 18 |
| TD023F | Binding region | 28 | 29 |
| TD031R | Binding region | 27 | 26 |
| TD006R | Full length | 40 | 39 |
| TD014R | Full length | 47 | 46 |

The three over-counts are also the three rows whose `GC%` is off by more than a
point, so one slip propagated. Anything measuring against this file must re-parse
the sequence rather than trust these cells.

`TD014R` and `TD031R` are also the only two rows whose Tm is more than 1.5 C from
what primer3 computes. Their Tm was probably read off the mis-parsed region. The
offset in `lib/primers.py` was fitted without those four rows.

## The `unsure` column is data, not a gate

`data/ytk_overhangs.tsv` still has an `unsure` column, and it is empty in every
row. Types 8 and 8a were the only ones ever marked, and the paper plus pYTK083-085
and pYTK089-091 confirmed both.

`ytk-add-overhangs` used to read that column and refuse to run, behind an
`--approve-unsure` flag. The flag was deleted in v0.5.0: with nothing in the
column it could never fire, and a flag that cannot fire is worse than no flag,
because the SKILL.md described a safeguard that was not there.

**So filling that cell in again will not stop anything.** The column is kept
because `tests/test_parts.py` checks it stays empty, which is what keeps the
verification on record. If a part type ever does need a gate, the code for it
has to come back too.

## Removing a cut site: one base is usually enough

`lib/silent.py` swaps a codon for a synonym to remove an internal BsmBI or BsaI
site. It is one of the two pieces of code here that change a sequence - `lib/codons.py`
is the other - and only `ytk-remove-cut-sites` calls it.

On both real genes in `tests/data/` it needs **one base**, at the third position
of a codon: `GAG -> GAA` at base 339 of Pi_fim_NCS_c1, and `GGT -> GGA` at base
411 of Pi_fim_NCS_c3. That is the usual case, because the genetic code is
redundant at the third base.

Two things it will not do, both deliberate:

- **It never touches the start or stop codon.** `TAA` and `TAG` both stop, so a
  naive synonym swap would happily change one for the other.
- **It requires each swap to reduce the total site count**, so a change cannot
  remove one site and create another somewhere else.

If no synonym removes a site, it says so rather than returning a sequence that
still holds one.

## Two published plasmids are cloned backwards

In **pYTK047** (the 234r GFP dropout) and **pYTK096** (the pre-assembled URA3
integration vector) the part sits in the reverse orientation, so the fragment
BsaI cuts out spans the numbering origin. `data/build_parts_table.py` handles
this. Anything new that cuts parts out of maps has to as well.

## Internal BsmBI sites in the connectors are deliberate

**pYTK002-008**, **pYTK050**, **pYTK067-073** and **pYTK096** each carry a BsmBI
site inside the part. Those are the level-2 multigene assembly sites, not
contamination. Do not flag them.

Every other part plasmid carries exactly one BsaI site per strand and no internal
BsmBI, so a part built with these flanks comes out with a clean single BsaI pair.

## The map has to start somewhere, and where matters

A circle has no first base. `ytk-clone` starts the map at base 1 of the backbone,
so the gene is never split across the start of the file and the same fragment
always gives the same map.

`data/ytk_parts.tsv` stores pYTK001 already turned to a point 40 bp inside ColE1,
so generated maps line up with the ones made by hand in SnapGene. A backbone read
straight from the kit starts at base 1 of the BsmBI overhang instead, which is
**inside the piece that drops out** and therefore gone from the finished plasmid.
When that happens the map starts at the first base of the piece that was kept.

**Do not "simplify" pYTK001 to just another file lookup.** `clone.load_backbone`
keeps an explicit first branch for it, and the measurements are why. The stored
copy and `reference/ytk_plasmids/pYTK001.gb` are the same 2676 bp circle, turned
to different first bases: base 1 of the reference file sits at position **177** of
the stored copy, so the stored copy is the reference turned by **2499** bases.
Base 1 of the reference file is the first base of the 1034 bp BsmBI dropout, so
the first 40 bases of that file are simply **not in** the finished plasmid
(measured: the kept piece is 1650 bp and does not hold them). The 40-base anchor
misses, the fallback fires, and every
entry-reaction map comes out rotated. `tests/test_ytk.py::TestPTP416` compares
byte for byte against a hand-drawn SnapGene map and would fail. With the branch,
`--backbone pYTK001` and no flag give byte-identical files, and a test pins that.
The other 95 published plasmids need no special handling.

## pydna counts sticky ends twice

A cut piece includes its 4-base overhangs, so two pieces joined measure eight
bases more than the circle they make. The arithmetic that holds:

    finished plasmid = backbone piece + insert - 8

For pYTK001: the vector is 2676 bp, the dropout contributes 1030 and the kept
backbone 1646, but `cut_backbone` measures 1650.

## Tm: why there is a constant in the code

`lib/primers.py` adds **1.31 C** to what primer3 returns. That is not biology and
not a fudge for the pair-balancing rule. It is the gap between two calculators of
the same number: primer3 offline, and the NEB website's Phusion calculator, which
is what the oligo sheet was filled in from. primer3 reads consistently lower.

Fitted once, over the 45 sound rows of `oligo_stocks`. After it, every one of
those 45 lands within 0.51 C of the NEB value. `tests/test_primers.py` fails if it
drifts. **Do not tune it per primer.**

## Codon optimisation is the weak part, and it is meant to look weak

`lib/codons.py` rewrites every codon as the most-used codon for its amino acid.
That is all it does. It is worth knowing why it does so little.

**A synthesis company's optimiser is better.** Twist, IDT and GenScript all
optimise at order time against their own synthesis-feasibility model, and they
will not build a sequence their own screen rejects. Argmax recoding is a worse
algorithm: it flattens the rare-codon structure real genes use for the initiation
ramp and for co-translational folding pauses, and there are published cases of
"optimised" genes expressing below native. `ytk-codon-optimise` exists for the
cases where you need the sequence before the order goes in, or need exactly which
codons moved on the record. Its SKILL.md says so in the opening paragraph, on
purpose.

**GC falls a long way, and nothing fixes it.** Yeast's preferred codons are
AT-rich. On the four real genes in `tests/data/` GC goes 47.3 -> 30.3, 51.2 ->
32.7, 48.3 -> 31.8, 53.8 -> 33.4. CAI goes 0.63 -> 1.00, 0.61 -> 1.00, 0.63 ->
1.00, 0.58 -> 0.99. Both numbers are printed side by side because the trade is
the whole story. Those figures are here so drift is visible if the table changes.

**The rewrite really can create a cut site, so step 1c is not optional.** Six
three-codon windows of yeast top codons hold one, all of them BsaI, and four of
those start with Trp: Trp-Ser-Gln recodes to `TGGTCTCAA`, which is `GGTCTC` plus a
base. `ytk-codon-optimise` exits 1 when it has done that, and
`test_the_rewrite_can_create_a_cut_site` keeps the case on record rather than
leaving it as a worry nobody checked. It is one base to clear afterwards.

**Repeats get worse, not better.** Identical peptide motifs recode to identical
DNA. One case is handled - the same codon is never used more than twice running,
so a poly-Q tract comes out `CAACAACAG` and not `CAACAACAA` - and a repeated motif
is not, because avoiding that needs a search over the whole sequence. Measured and
reported, never fixed.

## Why there is no DNA Chisel, yet

DNA Chisel would give GC targeting and repeat avoidance for free, and it was the
obvious choice until the detail. **Its randomness lives in numpy's global RNG, and
reproducibility across Python sessions is an open bug** -
github.com/Edinburgh-Genome-Foundry/DnaChisel issue 13: seeding fixes it inside
one session and not between them, with no maintainer answer.

The randomness only engages once the constraints you want it for go on.
`CodonOptimize(method="use_best_codon")` alone is a direct substitution and
deterministic; add `AvoidPattern`, `EnforceGCContent` or `UniquifyAllKmers` and it
becomes a local search. So adopting it today buys the dependency now and the bug
at the exact moment the constraints go on, for DNA somebody pays to synthesise.

**The seam is there instead.** The rewrite is behind `codons.recode(sequence)`.
Adding GC targeting means adding `lib/chisel.py` as a second backend and an
`--engine` flag, and no SKILL.md changes. The flag is deliberately absent until
there is a second engine to pick. Do not re-litigate this without re-reading
issue 13 first.

## There is no seed, and there must not be one

Weighted sampling from the usage table would match natural codon distribution
better than argmax. It is rejected anyway: it makes the output depend on a seed
somebody has to carry, on CPython's `random` staying stable, and on the order
values are drawn in. "There is no seed to remember" is a stronger guarantee than
any distribution argument, for a tool whose output gets ordered from a vendor.

What holds the output still, in order of how easily each could be broken:

- **Ties break alphabetically by codon.** The shares in the table are rounded to
  two places, so ties are real - human arginine has `AGA` and `AGG` both at 0.21.
  Without the rule the winner follows dict insertion order.
- **The package version is in the `codon_opt_method` column**, because the table
  lives in the package. A bump that moves one codon would otherwise change a
  sequence with nothing in the output saying why.
- **One golden sequence is pinned as a literal** in `tests/test_codons.py`.
- The loop runs left to right, never backtracks, and reads no set.

## Biopython already has CAI, and almost has the rewrite

`Bio.SeqUtils.CodonAdaptationIndex` implements Sharp & Li (1987). It is used for
CAI, so the numbers are comparable with anyone else's. Two things to know:

**The stop codons are deliberately left out of the index.** `calculate()` skips a
stop codon only when the index has no weight for it. An index built the usual way,
from whole genes, does carry stop weights, so it scores the terminator as though a
gene could have chosen a different one - and the same gene then scores differently
depending on whether it ends `TAA` or `TGA`. Sharp & Li exclude it.
`metrics._index` therefore builds from the 61 sense codons, and
`test_top_codons_score_one` is the check that this is right: with the stops in, an
all-top-codon sequence scored 0.956 instead of 1.0.

**`optimize()` is not used, though it is the same algorithm.** It replaces the
stop codon with whichever stop is most used, where a gene here keeps its own, and
on a tie it either raises or picks silently rather than by a stated rule. So
`lib/codons.py` mirrors it in about fifteen lines, and
`test_it_agrees_with_biopython_on_the_interior` pins the two together on the
codons where they should agree. `optimize()` needs stop weights, so that test
builds its own full index rather than weakening the one CAI uses.

## CAI cannot tell you which organism a gene was optimised for

This looked like a free feature - measure CAI against yeast, E. coli and human,
and the highest one names the source. It is wrong, and `input.csv` carries one CAI
column for that reason.

CAI is a geometric mean of each codon's share of its own family, so **a table with
more even codon usage scores every sequence higher.** The human table is flatter
than the yeast one: mean weight 0.73 against 0.66. A random ORF that nothing
optimised scores 0.638 against human and 0.599 against yeast - the same ranking
the four real fungal genes in `tests/data/` get. The columns would have measured
the shape of the tables and looked like they measured the gene.

Comparing codon usage between organisms properly means comparing the usage
frequencies or RSCU directly, not CAI. Nothing here needs that, because the only
question this plugin has to answer is whether a gene is already yeast-adapted.
