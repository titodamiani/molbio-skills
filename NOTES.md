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
site. It is the only code here that changes a sequence, and only
`ytk-cds-qc --remove-sites` calls it.

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
