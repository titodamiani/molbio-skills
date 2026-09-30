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
import enzymes  # noqa: E402
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
FRAGMENTS = {name: flanks.flank(seq, flanks.adapters("3"))
             for name, seq in GENES.items()}


def build(gene_name):
    return clone.assemble(FRAGMENTS[gene_name], BACKBONE, clone.BsmBI, "3")


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
            clone.assemble(GENES["Pi_fim_NCS_c5"], BACKBONE, clone.BsmBI, "3")

    def test_the_wrong_part_type_is_refused(self):
        # Type 5 flanks cut cleanly, but the overhangs do not fit pYTK001.
        fragment = flanks.flank(GENES["Pi_fim_NCS_c5"], flanks.adapters("5"))
        with self.assertRaises(clone.WrongFragment):
            clone.assemble(fragment, BACKBONE, clone.BsmBI, "3")

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
                clone.assemble(fragment, BACKBONE, clone.BsmBI, "3")
        finally:
            clone.cut_insert = real_cut_insert

        # and the honest case still works once the meddling is undone
        self.assertIsNotNone(
            sg.find_in_circle(clone.assemble(fragment, BACKBONE, clone.BsmBI, "3"),
                              gene))

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
        self.assertTrue(enzymes.beyond_flanks(FRAGMENTS["Pi_fim_NCS_c1"]))

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
        expected = (len(clone.cut_backbone(BACKBONE, clone.BsmBI))
                    + len(clone.cut_insert(FRAGMENTS["Pi_fim_OMT_c1"], clone.BsmBI))
                    - 2 * 4)
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
        """One folder per backbone under plasmids/, so a run against another
        vector cannot be mistaken for an entry-vector run."""
        return self.out / "plasmids" / clone.BACKBONE_NAME

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
        self.assertEqual("plasmid,sequence,part_type,codon_opt,notes", rows[0])
        self.assertTrue(rows[1].startswith("pTP412,Pi_fim_NCS_c1,"))

    def test_the_format_follows_the_fragments_with_no_flag(self):
        """The format is decided once, when the fragments are written. Reading it
        back off those files is what stops a run coming out half GenBank and
        half .dna."""
        (Path(self.tmp.name) / "already_here.dna").write_bytes(b"")
        self.run_clone()
        self.assertEqual(["pTP412.dna", "summary.csv"],
                         sorted(p.name for p in self.maps().iterdir()))


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
        self.assertEqual("pTP412", clone.plasmid_file_name("my_gene", "pTP412", "pYTK001"))

    def test_the_fallback_name_says_which_vector(self):
        # A folder full of <gene>_plasmid.dna files does not say what they are
        # in. Naming the vector means the file still reads clearly later.
        self.assertEqual("my_gene_pYTK001",
                         clone.plasmid_file_name("my_gene", None, "pYTK001"))

    def test_the_fallback_names_the_backbone_actually_used(self):
        # The name would be a lie if it named the default vector after a run
        # against a different one.
        self.assertEqual(BACKBONE, clone.backbone_sequence())
        self.assertEqual("g_pYTK047",
                         clone.plasmid_file_name("g", None, "pYTK047"))


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




class TestTheNote(unittest.TestCase):
    """The notes column is built from facts and never read back apart.

    The old version wrote prose and then parsed it again, which needed a
    sentinel for "nothing to say" and a rule for carrying hand-edited clauses
    through. A hand edit that changed the wording could then be misread, and the
    all-clear was printed over a gene that really had a site taken out. There is
    nothing to misread now: the facts travel in their own columns.
    """

    def test_nothing_done_and_nothing_left_gives_an_empty_note(self):
        self.assertEqual("", notes.summarise(set(), set()))

    def test_what_was_done_and_what_is_left_both_appear(self):
        self.assertEqual("BsmBI site removed; BsaI site in the CDS",
                         notes.summarise({"BsmBI"}, {"BsaI"}))

    def test_codon_optimisation_names_its_method(self):
        self.assertEqual("codon_opt (JCat); BsmBI site removed",
                         notes.summarise({"BsmBI"}, set(), "JCat"))

    def test_the_removed_column_survives_a_round_trip(self):
        self.assertEqual({"BsaI", "BsmBI"},
                         notes.unpack(notes.pack({"BsmBI", "BsaI"})))

    def test_an_empty_removed_column_means_nothing_was_removed(self):
        self.assertEqual(set(), notes.unpack(""))


class TestWhatWasDoneComesFromTheInputFile(unittest.TestCase):
    """clone.py used to copy the notes column along the chain, and read it from
    the first --input file only, so a second fragment file's genes lost the
    record of what was cleared out of them. It now reads the facts from --genes,
    which is one file, so there is no second file to forget."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.out = root / "out"
        self.first = root / "first.csv"
        self.second = root / "second.csv"
        self.genes = root / "input.csv"
        self.first.write_text(
            "name,sequence\n"
            f"Pi_fim_NCS_c1,{FRAGMENTS['Pi_fim_NCS_c1']}\n")
        self.second.write_text(
            "name,sequence\n"
            f"Pi_fim_NCS_c3,{FRAGMENTS['Pi_fim_NCS_c3']}\n")
        self.genes.write_text(
            "name,plasmid,sequence,removed,codon_opt_method\n"
            f"Pi_fim_NCS_c1,pTP412,{GENES['Pi_fim_NCS_c1']},BsmBI,\n"
            f"Pi_fim_NCS_c3,pTP414,{GENES['Pi_fim_NCS_c3']},BsaI,\n")

    def tearDown(self):
        self.tmp.cleanup()

    def test_every_fragment_file_keeps_what_was_done_to_its_genes(self):
        argv = sys.argv
        sys.argv = ["clone.py", "--input", str(self.first), str(self.second),
                    "--genes", str(self.genes), "--outdir", str(self.out)]
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                clone.main()
        finally:
            sys.argv = argv

        summary = self.out / "plasmids" / clone.BACKBONE_NAME / "summary.csv"
        rows = {row["sequence"]: (row["plasmid"], row["notes"])
                for row in csv.DictReader(summary.open())}
        # Both of these genes still hold a site of their own. That half is
        # measured here rather than carried, so it shows up beside what was
        # removed instead of being hidden by a stale column.
        self.assertEqual(("pTP412", "BsmBI site removed; BsmBI site in the CDS"),
                         rows["Pi_fim_NCS_c1"])
        self.assertEqual(("pTP414", "BsaI site removed; BsaI site in the CDS"),
                         rows["Pi_fim_NCS_c3"])


class TestTheBackboneIsAnInput(unittest.TestCase):
    """--backbone takes a name or a file, and a whole vector is loaded as a
    whole vector."""

    def test_a_published_name_is_found_without_a_path(self):
        name, sequence = clone.load_backbone("pYTK047")
        self.assertEqual("pYTK047", name)
        self.assertGreater(len(sequence), 2000)

    def test_a_vector_is_never_read_as_one_of_its_features(self):
        # pYTK032 holds exactly one feature that looks like a gene, so the
        # feature picker used to hand back that gene, 714 bp, and call it the
        # backbone. A backbone is the whole plasmid.
        name, sequence = clone.load_backbone("reference/ytk_plasmids/pYTK032.gb")
        self.assertEqual("pYTK032", name)
        self.assertEqual(2382, len(sequence))

    def test_the_default_and_the_name_give_the_same_vector(self):
        self.assertEqual(clone.load_backbone(None), clone.load_backbone("pYTK001"))
        self.assertEqual(clone.backbone_sequence(), clone.load_backbone("pYTK001")[1])

    def test_an_unknown_name_stops_the_run(self):
        with self.assertRaises(SystemExit):
            clone.load_backbone("pYTK999")


class TestThePartTypeIsDeclared(unittest.TestCase):
    """--type is checked on every run, whatever --backbone and --enzyme say."""

    def fragment(self, gene, part_type):
        return flanks.flank(GENES[gene], flanks.adapters(part_type))

    def test_a_type_5_part_declared_as_type_3_is_refused(self):
        # This is the regression the flags existed to hide: with the check off,
        # a Type 5 fragment cloned into pYTK001 and wrote a map.
        with self.assertRaises(clone.WrongFragment):
            clone.assemble(self.fragment("Pi_fim_NCS_c5", "5"), BACKBONE,
                           clone.BsmBI, "3")

    def test_a_type_5_part_declared_as_type_5_builds(self):
        # The entry vector takes any part type by design. Declaring the truth
        # lets it through.
        plasmid = clone.assemble(self.fragment("Pi_fim_NCS_c5", "5"), BACKBONE,
                                 clone.BsmBI, "5")
        self.assertIsNotNone(sg.find_in_circle(plasmid, GENES["Pi_fim_NCS_c5"]))

    def test_the_junctions_read_the_same_through_either_enzyme(self):
        fragment = self.fragment("Pi_fim_NCS_c1", "3")
        outer = clone.insert_junctions(clone.cut_insert(fragment, clone.BsmBI),
                                       clone.BsmBI)
        inner = clone.insert_junctions(clone.cut_insert(fragment, clone.BsaI),
                                       clone.BsaI)
        self.assertEqual(outer, inner)
        self.assertEqual(flanks.junctions("3"), outer)

    def test_a_part_that_does_not_fit_the_backbone_says_so(self):
        # A real Type 3 part, so the type check passes. pYTK047 is a type 234
        # slot, so the two overhang pairs do not meet. That has to be a
        # sentence, not a traceback out of pydna.
        backbone = clone.load_backbone("pYTK047")[1]
        with self.assertRaises(clone.WrongFragment) as refused:
            clone.assemble(self.fragment("Pi_fim_NCS_c1", "3"), backbone,
                           clone.BsaI, "3")
        message = str(refused.exception)
        for overhang in ("TATG", "ATCC", "AACG", "GCTG"):
            self.assertIn(overhang, message)


class TestTheFlagsReachTheRun(unittest.TestCase):
    """The same checks, through main(), because that is where the flags used to
    switch them off."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name) / "out"

    def tearDown(self):
        self.tmp.cleanup()

    def fragments(self, part_type, gene="Pi_fim_NCS_c1", name="Pi_fim_NCS_c1"):
        path = Path(self.tmp.name) / f"fragments_{part_type}.csv"
        flanked = flanks.flank(GENES[gene], flanks.adapters(part_type))
        path.write_text("name,sequence,part_type\n"
                        f"{name},{flanked},{part_type}\n")
        return path

    def run_clone(self, fragments, *extra, outdir=None):
        argv = sys.argv
        sys.argv = ["clone.py", "--input", str(fragments),
                    "--outdir", str(outdir or self.out)] + list(extra)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                clone.main()
        finally:
            sys.argv = argv

    def test_naming_the_default_vector_changes_nothing(self):
        named = Path(self.tmp.name) / "named"
        fragments = self.fragments("3")
        self.run_clone(fragments)
        self.run_clone(fragments, "--backbone", "pYTK001", outdir=named)
        default = (self.out / "plasmids" / "pYTK001" / "Pi_fim_NCS_c1_pYTK001.gb")
        self.assertEqual(default.read_bytes(),
                         (named / "plasmids" / "pYTK001"
                          / "Pi_fim_NCS_c1_pYTK001.gb").read_bytes())

    def test_a_backbone_flag_does_not_switch_the_type_check_off(self):
        with self.assertRaises(SystemExit):
            self.run_clone(self.fragments("5"), "--type", "3",
                           "--backbone", "reference/ytk_plasmids/pYTK001.gb")
        self.assertFalse((self.out / "plasmids").exists())

    def test_a_fragment_file_that_disagrees_with_the_type_stops_the_run(self):
        with self.assertRaises(SystemExit):
            self.run_clone(self.fragments("4"), "--type", "3")
        self.assertFalse((self.out / "plasmids").exists())

    def test_the_folder_and_the_file_name_follow_the_backbone(self):
        # The same entry vector under another file name, so the reaction works
        # and only the name differs. pYTK047 cannot be used here: it accepts a
        # whole 2-3-4 cassette, which is more than one fragment.
        vector = Path(self.tmp.name) / "my_vector.gb"
        vector.write_bytes(Path("reference/ytk_plasmids/pYTK001.gb").read_bytes())
        self.run_clone(self.fragments("3"), "--backbone", str(vector))
        folder = self.out / "plasmids" / "my_vector"
        self.assertTrue(folder.is_dir())
        self.assertFalse((self.out / "plasmids" / "pYTK001").exists())
        self.assertIn("Pi_fim_NCS_c1_my_vector.gb",
                      [path.name for path in folder.iterdir()])


class TestThePartTypeTables(unittest.TestCase):
    """The two data tables have to agree, or a part type means one thing in the
    flanks and another in the check."""

    def test_every_type_has_the_published_overhangs_in_its_adapters(self):
        for part_type in flanks._table(flanks.PART_TYPE_TABLE):
            with self.subTest(part_type):
                adapters = flanks.adapters(part_type)
                upstream, downstream = flanks.junctions(part_type)
                left = (adapters["left_adapter"] + "ATG").upper()
                self.assertTrue(left.startswith(upstream))
                self.assertTrue(adapters["right_adapter"].upper().endswith(downstream))

    def test_a_missing_column_reads_as_an_empty_string(self):
        self.assertEqual("", flanks.adapters("1")["coding"])

    def test_an_unknown_type_stops_the_run(self):
        with self.assertRaises(SystemExit):
            flanks.adapters("nonsense")

    def test_a_type_with_no_published_pair_stops_the_run(self):
        # custom has NNNN adapters, so it has to be finished by hand.
        with self.assertRaises(SystemExit):
            flanks.junctions("custom")


class TestTheWholeFixingChain(unittest.TestCase):
    """ytk-cds-qc, then ytk-codon-optimise, then ytk-remove-cut-sites, then step 2.

    Closes the loop the empty codon_opt columns were left open for: the method
    string has to reach fragments/summary.csv, and the two fixing skills have to
    write files step 2 can read without knowing which of them ran.
    """

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.out = Path(self.tmp.name)
        self.step("ytk-cds-qc", "cds_qc.py", "--input", str(DATA / "genes.csv"))
        self.step("ytk-codon-optimise", "codon_optimise.py",
                 "--input", str(self.out / "input.csv"))
        self.step("ytk-remove-cut-sites", "remove_cut_sites.py",
                 "--input", str(self.out / "optimised_genes" / "input_optimised.csv"))
        self.step("ytk-add-overhangs", "add_overhangs.py", "--type", "3",
                 "--input", str(self.out / "corrected_genes" / "input_corrected.csv"))

    def tearDown(self):
        self.tmp.cleanup()

    def step(self, skill, script, *args):
        """Run one skill's script, quietly. Its exit code is not the point here.

        ytk-cds-qc exits 1 on a gene with a cut site, which two of these have,
        and that is the behaviour the chain exists to deal with rather than a
        failure of the chain.
        """
        module = load(ROOT / "skills" / skill / "scripts" / script)
        argv = sys.argv
        sys.argv = [script, "--outdir", str(self.out)] + list(args)
        try:
            with contextlib.redirect_stdout(io.StringIO()):
                module.main()
        except SystemExit:
            pass
        finally:
            sys.argv = argv

    def rows(self, *parts):
        with open(self.out.joinpath(*parts), newline="") as fh:
            return list(csv.DictReader(fh))

    def test_the_method_reaches_the_synthesis_order(self):
        for row in self.rows("fragments", "summary.csv"):
            with self.subTest(row["name"]):
                self.assertEqual("true", row["codon_opt"])
                self.assertTrue(row["notes"].startswith("codon_opt (argmax/4932/"))

    def test_a_site_the_rewrite_cleared_is_still_on_the_record(self):
        """It cannot be measured off the new sequence, so it has to be carried."""
        notes_by_gene = {row["name"]: row["notes"]
                         for row in self.rows("fragments", "summary.csv")}
        self.assertIn("BsmBI site removed", notes_by_gene["Pi_fim_NCS_c1"])
        self.assertIn("BsaI site removed", notes_by_gene["Pi_fim_NCS_c3"])

    def test_the_two_fixing_skills_write_the_same_shape(self):
        optimised = self.rows("optimised_genes", "input_optimised.csv")
        corrected = self.rows("corrected_genes", "input_corrected.csv")
        self.assertEqual(list(optimised[0]), list(corrected[0]))

    def test_the_plasmid_names_survive_both_fixes(self):
        plasmids = [row["plasmid"]
                    for row in self.rows("corrected_genes", "input_corrected.csv")]
        self.assertIn("pTP412", plasmids)

    def test_the_measurements_are_retaken_after_the_rewrite(self):
        before = {r["name"]: r for r in self.rows("input.csv")}
        after = {r["name"]: r
                 for r in self.rows("optimised_genes", "input_optimised.csv")}
        for name in before:
            with self.subTest(name):
                # Yeast's favourite codons are AT-rich, so GC always falls.
                self.assertLess(float(after[name]["gc"]), float(before[name]["gc"]))
                self.assertGreater(float(after[name]["cai_scer"]),
                                   float(before[name]["cai_scer"]))

    def test_every_changed_codon_is_on_the_record(self):
        """The console prints a table only for a small gene, so the file is the record."""
        changes = self.rows("optimised_genes", "changes.csv")
        self.assertTrue(changes)
        for row in changes:
            with self.subTest(row["name"]):
                self.assertNotEqual(row["was"], row["now"])

    def test_nothing_still_holds_a_cut_site(self):
        for row in self.rows("corrected_genes", "input_corrected.csv"):
            with self.subTest(row["name"]):
                self.assertEqual([], enzymes.in_sequence(row["sequence"]))


if __name__ == "__main__":
    unittest.main(verbosity=2)
