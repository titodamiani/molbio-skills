"""Tests for reading the four input formats, and for asking once in a batch.

    python3 tests/test_inputs.py
"""
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "tests" / "data"
PLASMIDS = ROOT / "reference" / "ytk_plasmids"
sys.path.insert(0, str(ROOT / "lib"))
import cds          # noqa: E402
import report       # noqa: E402
import sequences    # noqa: E402


class TestTheFourFormats(unittest.TestCase):
    def test_fasta_and_csv_agree(self):
        fasta = [(name, seq) for name, seq, _ in sequences.read(DATA / "genes.fasta")]
        csv_rows = [(name, seq) for name, seq, _ in sequences.read(DATA / "genes.csv")]
        self.assertEqual(fasta, csv_rows)

    def test_only_the_csv_carries_plasmid_names(self):
        self.assertEqual([p for _, _, p in sequences.read(DATA / "genes.fasta")],
                         [None] * 4)
        self.assertIn("pTP412",
                      [p for _, _, p in sequences.read(DATA / "genes.csv")])

    def test_a_snapgene_map_gives_the_same_gene_as_the_fasta(self):
        genes = {name: seq for name, seq, _ in sequences.read(DATA / "genes.fasta")}
        [(name, seq, _)] = sequences.read(DATA / "pTP412.dna")
        self.assertEqual(seq, genes[name])

    def test_genbank_is_read(self):
        [(name, seq, _)] = sequences.read(PLASMIDS / "pYTK032.gb")
        self.assertEqual(name, "mTurquoise2")
        self.assertEqual(len(seq), 714)

    def test_an_unknown_format_is_refused(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "genes.txt"
            path.write_text("ATGAAATAA\n")
            with self.assertRaises(SystemExit):
                sequences.read(path)


class TestPickingTheCodingSequence(unittest.TestCase):
    def test_a_named_feature_wins(self):
        [(name, seq, _)] = sequences.read(DATA / "pTP412.dna", feature="CamR")
        self.assertEqual(name, "CamR")
        self.assertEqual(len(seq), 660)

    def test_a_name_that_is_not_there_says_what_is(self):
        with self.assertRaises(SystemExit) as stopped:
            sequences.read(DATA / "pTP412.dna", feature="nonsense")
        self.assertIn("CamR", str(stopped.exception))

    def test_the_marker_is_set_aside_not_returned(self):
        [(name, _, _)] = sequences.read(DATA / "pTP414.dna")
        self.assertEqual(name, "Pi_fim_NCS_c3")

    def test_several_candidates_stop_and_list_them(self):
        """pYTK096 is a pre-assembled vector, not a part plasmid, so several
        features survive the marker filter and it cannot be guessed."""
        with self.assertRaises(sequences.AmbiguousCDS) as ambiguous:
            sequences.read(PLASMIDS / "pYTK096.gb")
        table = ambiguous.exception.table()
        self.assertIn("name", table)
        self.assertGreater(len(ambiguous.exception.candidates), 1)

    def test_the_longest_candidate_is_offered_first(self):
        with self.assertRaises(sequences.AmbiguousCDS) as ambiguous:
            sequences.read(PLASMIDS / "pYTK096.gb")
        lengths = [f["end"] - f["start"] for f in ambiguous.exception.candidates]
        self.assertEqual(lengths, sorted(lengths, reverse=True))


class TestTheCodingSequenceCheck(unittest.TestCase):
    def test_a_good_cds_has_no_problems(self):
        self.assertEqual(cds.problems("ATGAAACGTTAA"), [])

    def test_every_fault_is_reported_not_just_the_first(self):
        problems = cds.problems("CCCAAACGTAAAA")
        self.assertEqual(len(problems), 3)
        self.assertTrue(any("ATG" in p for p in problems))
        self.assertTrue(any("stop codon" in p for p in problems))
        self.assertTrue(any("codons" in p for p in problems))

    def test_odd_letters_are_named(self):
        self.assertIn("N", cds.problems("ATGNNNTAA")[0])

    def test_a_stop_in_the_middle_is_a_warning_not_a_fault(self):
        sequence = "ATGTAACGTTAA"
        self.assertEqual(cds.problems(sequence), [])
        self.assertTrue(any("before the end" in w for w in cds.warnings(sequence)))

    def test_a_cut_site_is_a_warning(self):
        warnings = cds.warnings("ATG" + "CGTCTC" + "AAATAA")
        self.assertTrue(any("BsmBI" in w for w in warnings))


class TestAskingOnceForAWholeBatch(unittest.TestCase):
    def setUp(self):
        self.log = report.Report()

    def test_nothing_wrong_means_nothing_to_say(self):
        self.assertEqual(self.log.table(), "")
        self.assertEqual(self.log.question(), "")

    def test_one_table_holds_both_kinds(self):
        self.log.block("gene_a", "no stop codon")
        self.log.warn("gene_b", "holds a BsmBI site")
        table = self.log.table()
        self.assertIn("gene_a", table)
        self.assertIn("gene_b", table)
        self.assertIn("stops", table)
        self.assertIn("warning", table)

    def test_a_sequence_with_several_faults_is_named_once(self):
        self.log.block("gene_a", "no stop codon")
        self.log.block("gene_a", "not a whole number of codons")
        self.log.block("gene_b", "no ATG")
        self.assertEqual(self.log.blocked_names(), ["gene_a", "gene_b"])
        self.assertIn("2 of these", self.log.question())

    def test_warnings_alone_do_not_raise_a_question(self):
        self.log.warn("gene_a", "holds a BsaI site")
        self.assertNotEqual(self.log.table(), "")
        self.assertEqual(self.log.question(), "")


if __name__ == "__main__":
    unittest.main(verbosity=2)
