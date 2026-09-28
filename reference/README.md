# Reference material

Source files the code is checked against. Nothing here is imported at runtime
except by `data/build_parts_table.py`.

## `ytk_plasmids/`

The 96 published MoClo Yeast Toolkit plasmids, `pYTK001` to `pYTK096`, plus
`YTK_Parts.xls`, the published parts list.

- **Source:** Addgene kit #1000000061.
- **Paper:** Lee MW, DeLoache WC, Cervantes B, Dueber JE (2015). *A Highly
  Characterized Yeast Toolkit for Modular, Multipart Assembly.* ACS Synthetic
  Biology 4(9), 975-986. [doi:10.1021/sb500366v](https://doi.org/10.1021/sb500366v)

These are here so the repo checks itself. `data/ytk_parts.tsv` is generated from
them, and `tests/test_parts.py` rebuilds it and compares, so the table and the
generator cannot drift apart.

The two paper PDFs are **not** committed: 9 MB, and git keeps them forever. The
1 KB of the SI that the code actually needs is in `data/ytk_part_types.tsv`, with
the DOI in its header.

## `Plasmid_Generator.xlsx`

The spreadsheet the flanks were built in by hand before this plugin existed. Kept
because it is where the scaffold strings come from, and 12 KB is nothing.

It has two faults, both recorded in `../NOTES.md`:

1. The Type 3 and Type 3b right pieces are stored as `TAGATCC`, which adds a
   second stop codon to a CDS that already has one. `data/ytk_overhangs.tsv`
   stores plain `ATCC`.
2. Its `Left` and `Right` formulas skip two columns, so the advanced path
   silently drops adapters for types 2, 3a, 4 and 4a. Type 3 is unaffected.
