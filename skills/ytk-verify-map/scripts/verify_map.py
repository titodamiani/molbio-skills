"""Check SnapGene .dna files made by this plugin.

    python3 qc.py out/pTP416.dna --reference refs/pTP416.dna
    python3 qc.py out/*.dna --genes genes.fasta

Without a reference it checks the file on its own: the flags, the gene's
reading frame, and that the start of the map does not fall inside the gene.

With a reference it also compares the two as circles. A plasmid is a loop, so
two files can hold the same DNA while starting at different points. Comparing
them as plain text would wrongly report a difference.

It prints one line per check and ends with a count. Exit code 1 means
something failed.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "lib"))
import snapgene as sg

STOP_CODONS = ("TAA", "TAG", "TGA")


def check_file(path, reference=None, genes=()):
    """Run every check on one file. Returns a list of (passed, message)."""
    record = sg.read_dna(path)
    sequence = record["sequence"]
    results = [
        (record["double_stranded"], "double strand bit is set"),
        (len(sequence) > 0, f"length is {len(sequence)} bp"),
    ]

    if reference:
        wanted = sg.read_dna(reference)["sequence"]
        if record["circular"]:
            turn = sg.rotation_of(wanted, sequence)
            if turn is None:
                results.append((False, "does not match the reference"))
            elif turn == 0:
                results.append((True, "matches the reference exactly"))
            else:
                results.append((True, f"matches the reference, turned by {turn} bp"))
        else:
            results.append((sequence == wanted,
                            "matches the reference" if sequence == wanted
                            else "does not match the reference"))

    found = [(name, seq) for name, seq in genes
             if sg.find_in_circle(sequence, seq) is not None]

    if record["circular"]:
        results.append((True, "file is circular"))
        for name, gene in found:
            at = sg.find_in_circle(sequence, gene)
            results.append((at + len(gene) <= len(sequence),
                            f"{name} is not split by the start of the map"))
            results.append((len(gene) % 3 == 0,
                            f"{name} length {len(gene)} bp is a whole number of codons"))
            results.append((gene.startswith("ATG"), f"{name} starts with ATG"))
            results.append((gene[-3:] in STOP_CODONS,
                            f"{name} ends with a stop codon"))
        if genes and not found:
            results.append((False, "no gene from the gene file was found"))
    else:
        results.append((not record["circular"], "file is linear"))
        if genes:
            results.append((any(sequence == seq for _, seq in genes),
                            "sequence is unchanged from the input"))

    return results


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("files", nargs="+", help=".dna files to check")
    ap.add_argument("--reference", help="a hand-made .dna file to compare against")
    ap.add_argument("--genes", help="FASTA or CSV file of the cloned genes")
    args = ap.parse_args()

    if args.reference and len(args.files) > 1:
        sys.exit("--reference works with one file at a time")

    genes = [(n, s) for n, s, _ in sg.read_genes(args.genes)] if args.genes else []
    failed = 0

    for path in args.files:
        print(f"\n{Path(path).name}")
        for passed, message in check_file(path, args.reference, genes):
            print(f"  {'ok  ' if passed else 'FAIL'}  {message}")
            failed += not passed

    print(f"\n{failed} failed" if failed else "\nall checks passed")
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
