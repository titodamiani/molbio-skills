"""Where a run writes its results.

One rule, shared by every skill that writes anything, so a batch run cannot
scatter its parts across four different places:

    the folder asked for, or a ytk_output/ folder beside the input file

Output is always a folder and never a loose file, whatever the input was. A
single pasted sequence gets the same tree as a file of thirty genes, so there is
only ever one shape to look for afterwards.

There is no case here for an input that is not a file on disk. A skill that is
handed sequences in the chat writes them to a CSV first and then runs the normal
path, which keeps one code path rather than two.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import metrics

FOLDER_NAME = "ytk_output"
INPUT_NAME = "input.csv"


def folder(outdir, input_path):
    """The output folder, made if it is not there yet.

    Says so when it picked the folder itself. A default that is never printed is
    how someone ends up with results they cannot find.
    """
    if outdir:
        chosen = Path(outdir)
    else:
        chosen = _beside(Path(input_path).resolve())
        print(f"No output folder given, so everything goes in {chosen}")
    chosen.mkdir(parents=True, exist_ok=True)
    return chosen


def _beside(input_path):
    """The folder to use for an input that named none.

    Steps run by hand feed each other: ytk-clone reads the summary ytk-add-
    overhangs wrote, which already sits inside the output folder. Putting a new
    ytk_output/ beside that file would bury one folder in the last one on every
    step, so an input that is already under one goes back into it.
    """
    for parent in input_path.parents:
        if parent.name == FOLDER_NAME:
            return parent
    return input_path.parent / FOLDER_NAME


# What ytk-remove-cut-sites and ytk-codon-optimise fill in. Written empty by
# ytk-cds-qc so input.csv, input_corrected.csv and input_optimised.csv are one
# shape, and a later step reads the same columns whichever of the three it got.
FIXED_BY_LATER_STEPS = ["removed", "codon_opt_method"]

HEADER = ["name", "plasmid", "sequence"] + metrics.COLUMNS + FIXED_BY_LATER_STEPS

SUMMARY_HEADER = ["name", "input_sequence", "new_sequence", "codon_opt", "notes"]


def _row(name, plasmid, sequence, removed="", codon_opt_method=""):
    """One gene row, measured from the sequence in that row.

    Measured here rather than passed in, so a file written after a fix carries
    the numbers for the fixed sequence and never the old ones. An empty cell
    means the number could not be taken, not that it was zero.
    """
    measured = metrics.measure(sequence)
    return ([name, plasmid or "", sequence]
            + ["" if measured[column] is None else measured[column]
               for column in metrics.COLUMNS]
            + [removed, codon_opt_method])


def _write(path, header, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        writer.writerows(rows)
    return path


def write_input_csv(rows, folder):
    """A copy of the input, measured, in the one shape every later step reads.

    `rows` is what sequences.read returns: (name, sequence, plasmid or None).
    Writing it here rather than copying the input file means a FASTA, a GenBank
    map or a .dna file reaches the rest of the run as a CSV with a plasmid
    column, so there is one input format downstream instead of four.

    The measurements ride along on the sequence's own row. Measured once, they
    tell whoever reads the file which fixing skill a gene needs, and a later
    note is rebuilt from the columns rather than from re-measuring or from
    reading a sentence back apart. Everything downstream picks columns out by
    header name, so the extra ones cost it nothing.
    """
    return _write(Path(folder) / INPUT_NAME, HEADER,
                  [_row(name, plasmid, sequence)
                   for name, sequence, plasmid in rows])


def write_genes_csv(path, rows):
    """The handoff a fixing skill writes, in the same shape as input.csv.

    `rows` is (name, plasmid, sequence, removed, codon_opt_method). Same shape as
    input.csv on purpose: whichever file the next step is pointed at, it reads
    the same columns, so nothing downstream has to know which fixes ran.
    """
    return _write(path, HEADER, [_row(*row) for row in rows])


def write_summary_csv(path, rows):
    """The file for reading: the old sequence beside the new one.

    `rows` is (name, input_sequence, new_sequence, codon_opt, notes). Separate
    from write_genes_csv because the two answer different questions - this one
    is read by a person, that one is read by the next step.
    """
    return _write(path, SUMMARY_HEADER, rows)
