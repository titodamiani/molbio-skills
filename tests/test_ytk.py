"""Tests for the ytk-skills plugin.

Run them like this:

    python3 tests/test_ytk.py

They need nothing installed beyond biopython and pydna, which the scripts
install by themselves. The reference .dna files in tests/data were made by
hand in SnapGene.

The five cases each cover something different:

    pTP412          the reference starts at a different point on the circle,
                    and the gene holds an extra BsmBI site inside it
    pTP414          the gene holds a BsaI site inside it, which must survive
    pTP416          an ordinary case, so the map should match exactly
    pTP0457         a long insert with no hand-made reference, so only the
                    properties of the map are checked
    Pi_fim_NCS_c1   the linear writer, which the plasmid cases never touch
"""
import contextlib
import csv
import importlib.util
import io
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DATA = Path(__file__).resolve().parent / "data"
sys.path.insert(0, str(ROOT / "lib"))
import flanks  # noqa: E402
import notes  # noqa: E402
import snapgene as sg  # noqa: E402


def same_circle(a, b):
    """True when two sequences describe the same circular DNA.

    A test helper, not shipped code: only the tests ever ask this question as a
    yes or no. Everything in lib/ wants the offset from sg.rotation_of instead.
    """
    return sg.rotation_of(a, b) is not None


def load(path):
    """Import a script from a folder whose name has a dash in it."""
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


clone = load(ROOT / "skills" / "ytk-clone" / "scripts" / "clone.py")
annotate = load(ROOT / "skills" / "ytk-annotate-map" / "scripts" / "annotate_map.py")
overhangs = load(ROOT / "skills" / "ytk-add-overhangs" / "scripts" / "add_overhangs.py")
verify = load(ROOT / "skills" / "ytk-verify-map" / "scripts" / "verify_map.py")

GENES = {name: seq for name, seq, _ in sg.read_genes(DATA / "genes.fasta")}
BACKBONE = clone.backbone_sequence()

# ytk-clone now takes the fragment as ordered, so the genes are flanked here
# the same way ytk-add-overhangs flanks them. That makes these cases a test of
# the two skills together.
ADAPTERS = overhangs.read_overhangs()
FRAGMENTS = {name: flanks.flank(seq, ADAPTERS["3"]) for name, seq in GENES.items()}


def build(gene_name):
    return clone.assemble(FRAGMENTS[gene_name], BACKBONE)


def reference(name):
    return sg.read_dna(DATA / f"{name}.dna")["sequence"]


class TestCircleCompare(unittest.TestCase):
    """A plasmid is a loop, so the same DNA can start at different points."""

    def test_a_turned_circle_is_still_the_same_circle(self):
        circle = "AAAACCCCGGGGTTTT"
        self.assertTrue(same_circle(circle, sg.rotate_to(circle, "GGGG")))

    def test_different_dna_is_not_the_same_circle(self):
        self.assertFalse(same_circle("AAAACCCC", "AAAACCCG"))

    def test_plain_text_compare_would_have_missed_it(self):
        turned = sg.rotate_to("AAAACCCCGGGGTTTT", "GGGG")
        self.assertNotEqual("AAAACCCCGGGGTTTT", turned)


class TestSequencesAreNeverChanged(unittest.TestCase):
    """The one promise that matters most: a gene goes in and comes out the
    same. Genes with a cut site inside them are cut and joined back on the way
    through, so this is where a change could slip in unseen."""

    def test_every_gene_is_in_its_plasmid_exactly_as_given(self):
        for name, gene in GENES.items():
            with self.subTest(gene=name):
                self.assertIsNotNone(sg.find_in_circle(build(name), gene))

    def test_the_linear_file_holds_the_gene_exactly_as_given(self):
        for name, gene in GENES.items():
            with self.subTest(gene=name):
                with tempfile.TemporaryDirectory() as tmp:
                    path = Path(tmp) / "gene.dna"
                    sg.write_dna(path, gene, circular=False)
                    self.assertEqual(gene, sg.read_dna(path)["sequence"])

    def test_a_bare_gene_is_refused(self):
        # The flanks are no longer added here, so a bare gene has nothing to
        # cut. Adding the flanks for the person would assume a fragment design
        # they may never have ordered.
        with self.assertRaises(clone.WrongFragment):
            clone.assemble(GENES["Pi_fim_NCS_c5"], BACKBONE)

    def test_the_wrong_part_type_is_refused(self):
        # Type 5 flanks cut cleanly, but the overhangs do not fit pYTK001.
        fragment = flanks.flank(GENES["Pi_fim_NCS_c5"], ADAPTERS["5"])
        with self.assertRaises(clone.WrongFragment):
            clone.assemble(fragment, BACKBONE)

    def test_the_run_stops_if_a_gene_came_out_changed(self):
        # Make the assembly quietly change one base in the middle of the gene.
        # The sticky ends come from the flanks, so the loop still closes and
        # the plasmid still looks perfectly normal. Only the check catches it.
        # Without the check this would be written to a .dna file and used.
        gene = GENES["Pi_fim_NCS_c5"]
        fragment = FRAGMENTS["Pi_fim_NCS_c5"]
        real_cut_insert = clone.cut_insert

        def cut_a_changed_gene(sequence, *enzyme):
            swapped = "C" if sequence[30] != "C" else "A"
            return real_cut_insert(sequence[:30] + swapped + sequence[31:],
                                   *enzyme)

        clone.cut_insert = cut_a_changed_gene
        try:
            with self.assertRaises(clone.SequenceChanged):
                clone.assemble(fragment, BACKBONE)
        finally:
            clone.cut_insert = real_cut_insert

        # and the honest case still works once the meddling is undone
        self.assertIsNotNone(
            sg.find_in_circle(clone.assemble(fragment, BACKBONE), gene))

    def test_a_sequence_with_odd_letters_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "bad.csv"
            path.write_text("gene1,ATGNNNTAA\n")
            with self.assertRaises(SystemExit):
                sg.read_genes(path)


class TestPTP412(unittest.TestCase):
    """The reference starts elsewhere on the circle, and the gene holds an
    extra BsmBI site inside it."""

    def setUp(self):
        self.built = build("Pi_fim_NCS_c1")
        self.wanted = reference("pTP412")

    def test_the_gene_really_does_hold_an_extra_site(self):
        self.assertTrue(clone.has_internal_site(FRAGMENTS["Pi_fim_NCS_c1"]))

    def test_same_circle_as_the_reference(self):
        self.assertTrue(same_circle(self.wanted, self.built))

    def test_the_reference_starts_somewhere_else(self):
        self.assertNotEqual(0, sg.rotation_of(self.wanted, self.built))

    def test_the_gene_came_through_whole(self):
        self.assertIsNotNone(
            sg.find_in_circle(self.built, GENES["Pi_fim_NCS_c1"]))


class TestPTP414(unittest.TestCase):
    """The gene holds a BsaI site inside it, which must be left alone."""

    def test_the_gene_really_does_hold_an_internal_site(self):
        self.assertIn("GGTCTC", GENES["Pi_fim_NCS_c3"])

    def test_same_circle_as_the_reference(self):
        self.assertTrue(same_circle(reference("pTP414"), build("Pi_fim_NCS_c3")))

    def test_the_gene_is_untouched(self):
        self.assertIsNotNone(
            sg.find_in_circle(build("Pi_fim_NCS_c3"), GENES["Pi_fim_NCS_c3"]))


class TestPTP416(unittest.TestCase):
    """An ordinary case. The map should match the reference exactly."""

    def test_matches_exactly(self):
        self.assertEqual(reference("pTP416"), build("Pi_fim_NCS_c5"))


class TestPTP0457(unittest.TestCase):
    """A long insert with no hand-made reference, so check properties only."""

    def setUp(self):
        self.gene = GENES["Pi_fim_OMT_c1"]
        self.built = build("Pi_fim_OMT_c1")

    def test_length_is_the_backbone_piece_plus_the_cut_insert(self):
        # The two pieces join at two sticky ends. Each end is 4 bases that the
        # two pieces share, so the loop is 8 bases shorter than the sum.
        expected = (len(clone.cut_backbone(BACKBONE))
                    + len(clone.cut_insert(FRAGMENTS["Pi_fim_OMT_c1"])) - 2 * 4)
        self.assertEqual(expected, len(self.built))

    def test_still_matches_the_map_the_original_code_made(self):
        # tests/data/pTP0457.dna was not drawn by hand. It is what the working
        # code produced, kept here so a change in the logic gets noticed.
        self.assertTrue(same_circle(reference("pTP0457"), self.built))

    def test_the_gene_is_a_whole_number_of_codons(self):
        self.assertEqual(0, len(self.gene) % 3)
        self.assertTrue(self.gene.startswith("ATG"))
        self.assertIn(self.gene[-3:], ("TAA", "TAG", "TGA"))

    def test_the_start_of_the_map_does_not_split_the_gene(self):
        at = sg.find_in_circle(self.built, self.gene)
        self.assertLessEqual(at + len(self.gene), len(self.built))


class TestLinearWriter(unittest.TestCase):
    """The plasmid cases never exercise the linear writer."""

    def setUp(self):
        self.gene = GENES["Pi_fim_NCS_c1"]
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "gene.dna"
        sg.write_dna(self.path, self.gene, circular=False)
        self.written = sg.read_dna(self.path)

    def tearDown(self):
        self.tmp.cleanup()

    def test_sequence_is_unchanged(self):
        self.assertEqual(self.gene, self.written["sequence"])

    def test_file_is_linear_not_circular(self):
        self.assertFalse(self.written["circular"])

    def test_double_strand_bit_is_set(self):
        self.assertTrue(self.written["double_stranded"])

    def test_matches_the_hand_made_reference(self):
        self.assertEqual(reference("Pi_fim_NCS_c1"), self.written["sequence"])


class TestGenBankRoundTrip(unittest.TestCase):
    """GenBank is the portable format, for colleagues without SnapGene. A map
    written that way has to hold exactly what the SnapGene writer would hold,
    or the two formats quietly disagree about the same plasmid."""

    def setUp(self):
        self.plasmid = build("Pi_fim_NCS_c1")
        self.tmp = tempfile.TemporaryDirectory()
        self.gb = Path(self.tmp.name) / "pTP412.gb"
        self.dna = Path(self.tmp.name) / "pTP412.dna"

    def tearDown(self):
        self.tmp.cleanup()

    def test_genbank_round_trips_the_sequence(self):
        sg.write_map(self.gb, self.plasmid, circular=True)
        written = sg.read_map(self.gb)
        self.assertEqual(self.plasmid, written["sequence"])
        self.assertTrue(written["circular"])

    def test_genbank_and_dna_hold_the_same_sequence(self):
        sg.write_map(self.gb, self.plasmid, circular=True)
        sg.write_map(self.dna, self.plasmid, circular=True)
        self.assertEqual(sg.read_map(self.dna)["sequence"],
                         sg.read_map(self.gb)["sequence"])

    def test_a_linear_genbank_file_is_not_circular(self):
        path = Path(self.tmp.name) / "gene.gb"
        sg.write_map(path, GENES["Pi_fim_NCS_c1"], circular=False)
        self.assertFalse(sg.read_map(path)["circular"])


class TestGenBankFeatures(unittest.TestCase):
    """A feature running past the end of the circle is the one most likely to
    come out wrong, because GenBank writes it as a join() of two pieces while
    SnapGene writes an end number lower than its start."""

    def setUp(self):
        self.plasmid = build("Pi_fim_NCS_c1")
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "labelled.gb"

    def tearDown(self):
        self.tmp.cleanup()

    def features_back(self):
        from Bio import SeqIO
        return SeqIO.read(str(self.path), "genbank").features

    def test_genbank_feature_labels_round_trip(self):
        sg.write_map(self.path, self.plasmid, circular=True, features=[
            {"name": "ColE1", "type": "rep_origin", "start": 10, "end": 110,
             "strand": 1, "wrap_end": 0}])
        back = self.features_back()
        self.assertEqual(["ColE1"], [f.qualifiers["label"][0] for f in back])
        self.assertEqual(["rep_origin"], [f.type for f in back])
        self.assertEqual(10, int(back[0].location.start))
        self.assertEqual(110, int(back[0].location.end))

    def test_a_feature_that_crosses_the_origin_survives_genbank(self):
        # Starts 40 bases before the end of the circle and runs 60 past it.
        start = len(self.plasmid) - 40
        sg.write_map(self.path, self.plasmid, circular=True, features=[
            {"name": "wraps", "type": "misc_feature", "start": start,
             "end": start + 100, "strand": 1, "wrap_end": 60}])
        covered = sorted({int(base) for base in self.features_back()[0].location})
        self.assertEqual(100, len(covered))
        for base in (start, len(self.plasmid) - 1, 0, 59):
            self.assertIn(base, covered)


class TestGenBankEndToEnd(unittest.TestCase):
    """ytk-clone has to write GenBank on request, and ytk-verify-map has to read
    it back, or step 6 of the workflow stops working for anyone who picks it."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"
        fragments = Path(self.tmp.name) / "fragments.csv"
        fragments.write_text(
            "name,sequence,plasmid\n"
            f"Pi_fim_NCS_c1,{FRAGMENTS['Pi_fim_NCS_c1']},pTP412\n")
        self.argv = ["clone.py", "--input", str(fragments),
                     "--outdir", str(self.out)]

    def tearDown(self):
        self.tmp.cleanup()

    def run_clone(self, *extra):
        argv = sys.argv
        sys.argv = self.argv + list(extra)
        try:
            # main() prints its own report, which would bury the test results.
            with contextlib.redirect_stdout(io.StringIO()):
                clone.main()
        finally:
            sys.argv = argv

    def maps(self):
        """The map folder is named after the backbone actually used, so a run
        against another vector cannot be mistaken for an entry-vector run."""
        return self.out / f"{clone.BACKBONE_NAME}_maps"

    def test_clone_writes_genbank_by_default(self):
        self.run_clone()
        self.assertEqual(["pTP412.gb", "summary.csv"],
                         sorted(p.name for p in self.maps().iterdir()))

    def test_clone_writes_snapgene_when_asked(self):
        self.run_clone("--format", "dna")
        self.assertEqual(["pTP412.dna", "summary.csv"],
                         sorted(p.name for p in self.maps().iterdir()))

    def test_both_formats_give_the_same_plasmid(self):
        self.run_clone()
        self.run_clone("--format", "dna")
        self.assertEqual(sg.read_map(self.maps() / "pTP412.dna")["sequence"],
                         sg.read_map(self.maps() / "pTP412.gb")["sequence"])

    def test_verify_map_reads_a_genbank_file(self):
        self.run_clone()
        results = verify.check_file(self.maps() / "pTP412.gb",
                                    genes=[("Pi_fim_NCS_c1", GENES["Pi_fim_NCS_c1"])])
        failed = [message for passed, message in results if not passed]
        self.assertEqual([], failed)

    def test_the_summary_carries_the_plasmid_name(self):
        self.run_clone()
        rows = (self.maps() / "summary.csv").read_text().splitlines()
        self.assertEqual("name,part_type,plasmid,notes", rows[0])
        self.assertIn("pTP412", rows[1])


class TestPartsTable(unittest.TestCase):
    """The parts table is generated. Two rows sharing a sequence means the
    generator picked the wrong fragment out of a plasmid, which is how three
    different bacterial markers once ended up as the same GFP dropout."""

    def test_no_two_rows_share_a_sequence(self):
        import collections
        import csv
        with open(ROOT / "data" / "ytk_parts.tsv", newline="") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        seen = collections.Counter(row["sequence"] for row in rows)
        shared = [seq for seq, count in seen.items() if count > 1]
        names = [row["name"] for row in rows if row["sequence"] in shared]
        self.assertEqual([], names)


class TestInput(unittest.TestCase):
    """FASTA and CSV must give the same answer."""

    def test_csv_and_fasta_hold_the_same_genes(self):
        self.assertEqual([(n, s) for n, s, _ in sg.read_genes(DATA / "genes.fasta")],
                         [(n, s) for n, s, _ in sg.read_genes(DATA / "genes.csv")])

    def test_a_csv_without_a_header_row_works_too(self):
        self.assertEqual([("gene1", "ATGAAATAA", None)],
                         self.read("gene1,ATGAAATAA\n"))

    def test_headers_are_matched_whatever_the_capitals(self):
        self.assertEqual([("gene1", "ATGAAATAA", None)],
                         self.read("NAME,SEQUENCE\ngene1,ATGAAATAA\n"))

    def test_a_byte_order_mark_from_excel_is_ignored(self):
        # A CSV saved by Excel starts with a BOM. Read as plain UTF-8 the first
        # header becomes "﻿name", matches nothing, and the header row is
        # then parsed as a gene.
        self.assertEqual([("gene1", "ATGAAATAA", None)],
                         self.read("﻿name,sequence\ngene1,ATGAAATAA\n"))

    def read(self, text, suffix=".csv"):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / f"genes{suffix}"
            path.write_text(text)
            return sg.read_genes(path)


class TestPlasmidNames(unittest.TestCase):
    """Real plasmid numbers are not a tidy series, so the names come from the
    gene file rather than being counted out by the script."""

    def read(self, text):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.csv"
            path.write_text(text)
            return sg.read_genes(path)

    def test_the_plasmid_column_is_used(self):
        self.assertEqual(
            [("a", "ATGAAATAA", "pTP412"), ("b", "ATGCCCTAA", "pTP768")],
            self.read("name,sequence,plasmid\na,ATGAAATAA,pTP412\nb,ATGCCCTAA,pTP768\n"))

    def test_the_numbers_do_not_have_to_be_in_a_series(self):
        genes = self.read(
            "name,sequence,plasmid\na,ATGAAATAA,pTP412\nb,ATGCCCTAA,pTP768\nc,ATGGGGTAA,pTP002\n")
        self.assertEqual(["pTP412", "pTP768", "pTP002"], [g[2] for g in genes])

    def test_a_headerless_csv_can_carry_a_plasmid_name_too(self):
        self.assertEqual([("a", "ATGAAATAA", "pTP412")],
                         self.read("a,ATGAAATAA,pTP412\n"))

    def test_a_blank_plasmid_cell_falls_back(self):
        self.assertIsNone(self.read("name,sequence,plasmid\na,ATGAAATAA,\n")[0][2])

    def test_fasta_never_carries_a_plasmid_name(self):
        for _, _, plasmid in sg.read_genes(DATA / "genes.fasta"):
            self.assertIsNone(plasmid)

    def test_the_real_plasmid_names_come_through(self):
        names = {n: p for n, _, p in sg.read_genes(DATA / "genes.csv")}
        self.assertEqual("pTP412", names["Pi_fim_NCS_c1"])
        self.assertEqual("pTP0457", names["Pi_fim_OMT_c1"])

    def test_a_repeated_plasmid_name_stops_the_run(self):
        # Without this, the second file would silently overwrite the first.
        with self.assertRaises(SystemExit):
            self.read("name,sequence,plasmid\na,ATGAAATAA,pTP412\nb,ATGCCCTAA,pTP412\n")

    def test_a_repeated_gene_name_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.read("name,sequence\na,ATGAAATAA\na,ATGCCCTAA\n")

    def test_a_given_name_is_used_as_it_stands(self):
        self.assertEqual("pTP412", clone.plasmid_file_name("my_gene", "pTP412"))

    def test_the_fallback_name_says_which_vector(self):
        # A folder full of <gene>_plasmid.dna files does not say what they are
        # in. Naming the vector means the file still reads clearly later.
        self.assertEqual("my_gene_pYTK001",
                         clone.plasmid_file_name("my_gene", None))

    def test_the_fallback_vector_is_the_one_actually_used(self):
        # The name would be a lie if it drifted from the backbone in the table.
        self.assertEqual(BACKBONE, clone.backbone_sequence())
        self.assertIn(clone.BACKBONE_NAME, clone.plasmid_file_name("g", None))


class TestAnnotation(unittest.TestCase):
    """Labels come from the shared parts table, not from another map file."""

    def setUp(self):
        self.built = build("Pi_fim_NCS_c5")
        self.features = annotate.build_features(
            self.built, annotate.load_parts(),
            [("Pi_fim_NCS_c5", GENES["Pi_fim_NCS_c5"])])
        self.named = {f["name"]: f for f in self.features}

    def test_the_backbone_parts_are_labelled(self):
        for name in ("ColE1", "CamR", "CamR Promoter", "CamR Terminator"):
            self.assertIn(name, self.named)

    def test_the_gene_is_labelled(self):
        self.assertIn("Pi_fim_NCS_c5", self.named)

    def test_cole1_is_764_bp_not_the_whole_plasmid(self):
        cole1 = self.named["ColE1"]
        self.assertEqual(764, cole1["end"] - cole1["start"] + cole1["wrap_end"])

    def test_every_label_sits_where_its_sequence_is(self):
        for feature in self.features:
            self.assertLess(feature["start"], len(self.built))


if __name__ == "__main__":
    unittest.main(verbosity=2)


class TestOrderNote(unittest.TestCase):
    """The notes column of fragments/summary.csv is what a synthesis order is
    placed from, so it must never claim a gene was clean when it was changed.

    The cell is a CSV field someone can edit in Excel. Matching it against a
    fixed phrase used to drop `BsmBI site removed.` for its full stop, and the
    all-clear was then printed over a gene that really had a site taken out.
    """

    CLEAN = "ATGAAACCCGGGTAA"
    DIRTY = "ATGCGTCTCAAATAA"      # holds a BsmBI site

    def test_a_hand_edited_note_is_never_turned_into_the_all_clear(self):
        for edited in ("BsmBI site removed.",
                       "BsmBI Site Removed",
                       "BsmBI site removed (checked by hand)"):
            with self.subTest(edited):
                self.assertNotEqual(notes.NO_SITES,
                                    notes.for_order(self.CLEAN, edited))

    def test_a_hand_edited_note_is_carried_through_word_for_word(self):
        self.assertEqual("BsmBI site removed.",
                         notes.for_order(self.CLEAN, "BsmBI site removed."))

    def test_only_an_empty_cell_gives_the_all_clear(self):
        self.assertEqual(notes.NO_SITES, notes.for_order(self.CLEAN, ""))
        self.assertEqual(notes.NO_SITES,
                         notes.for_order(self.CLEAN, notes.NO_SITES))

    def test_a_site_in_the_cds_is_measured_and_not_inherited(self):
        # Inheriting this half is what would let a stale note bless an unsafe
        # fragment, so it is always read off the sequence being ordered.
        self.assertEqual(f"BsmBI {notes.PRESENT}",
                         notes.for_order(self.DIRTY, ""))
        self.assertEqual(notes.NO_SITES,
                         notes.for_order(self.CLEAN, f"BsaI {notes.PRESENT}"))


class TestNotesAcrossSeveralInputFiles(unittest.TestCase):
    """clone.py used to read its notes column from the first --input file only,
    so a second file's genes lost the record of what was cleared out of them."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.out = root / "out"
        self.first = root / "first.csv"
        self.second = root / "second.csv"
        self.first.write_text(
            "name,sequence,plasmid,notes\n"
            f"Pi_fim_NCS_c1,{FRAGMENTS['Pi_fim_NCS_c1']},pTP412,"
            "BsmBI site removed\n")
        self.second.write_text(
            "name,sequence,plasmid,notes\n"
            f"Pi_fim_NCS_c3,{FRAGMENTS['Pi_fim_NCS_c3']},pTP414,"
            "BsaI site removed\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_every_input_file_keeps_its_notes(self):
        argv = sys.argv
        sys.argv = ["clone.py", "--input", str(self.first), str(self.second),
                    "--outdir", str(self.out)]
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                clone.main()
        finally:
            sys.argv = argv

        summary = (self.out / f"{clone.BACKBONE_NAME}_maps" / "summary.csv")
        rows = {row["name"]: row["notes"]
                for row in csv.DictReader(summary.open())}
        self.assertEqual({"Pi_fim_NCS_c1": "BsmBI site removed",
                          "Pi_fim_NCS_c3": "BsaI site removed"}, rows)
