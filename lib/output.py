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
from pathlib import Path

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


def write_input_csv(rows, folder):
    """A copy of the input, in the one shape every later step reads.

    `rows` is what sequences.read returns: (name, sequence, plasmid or None).
    Writing it here rather than copying the input file means a FASTA, a GenBank
    map or a .dna file reaches the rest of the run as a CSV with a plasmid
    column, so there is one input format downstream instead of four.
    """
    path = Path(folder) / INPUT_NAME
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["name", "plasmid", "sequence"])
        writer.writerows([name, plasmid or "", sequence]
                         for name, sequence, plasmid in rows)
    return path
