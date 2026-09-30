"""The codon usage table, the rewrite, and the measurements.

    python3 tests/test_codons.py

The table is the part most likely to be wrong, because it comes from a package
rather than from this repo, so it is tested hardest: a bump that swapped the
yeast table for another organism's would change every sequence this plugin
writes, and nothing else here would notice.
"""
import sys
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
import cds       # noqa: E402
import enzymes   # noqa: E402
import codons    # noqa: E402
import metrics   # noqa: E402
import silent    # noqa: E402
import snapgene as sg  # noqa: E402

import python_codon_tables as pct  # noqa: E402
from Bio.Seq import Seq  # noqa: E402
from Bio.SeqUtils import CodonAdaptationIndex  # noqa: E402

GENES = ROOT / "tests" / "data" / "genes.csv"

# A short gene, recoded, pinned exactly. Any change to the algorithm or to the
# table has to walk past this, which is the point: the sequences this writes get
# ordered from a synthesis company, so a silent change is the worst outcome.
GOLDEN_IN = "ATGGCTAGCAAGGGCGAAGAACTGTTCACCGGCGTTCAGCCGATTCTGGTGGAACTGGATTAA"
GOLDEN_OUT = "ATGGCTTCTAAAGGTGAAGAATTGTTTACTGGTGTTCAACCAATTTTGGTTGAATTGGATTAA"


def full_index(taxid=codons.YEAST):
    """A CAI index with the stop codons in it, which only optimize() needs.

    metrics._index leaves them out on purpose, so that calculate() skips the
    terminator as Sharp & Li do. optimize() reads the index the other way round
    and wants a weight for every codon, so the comparison builds its own.
    """
    index = CodonAdaptationIndex(["ATG"])
    index.clear()
    index.update({codon: share / max(family.values())
                  for family in pct.get_codons_table(taxid).values()
                  for codon, share in family.items()})
    return index


def gene(peptide, stop="TAA"):
    """A coding sequence for a peptide, using each amino acid's top codon."""
    best = codons._choices(codons.YEAST)
    return "ATG" + "".join(best[amino][0] for amino in peptide) + stop


class TestTheTable(unittest.TestCase):
    def test_every_codon_is_there_once(self):
        ranked = codons.ranked()
        listed = [codon for order in ranked.values() for codon in order]
        self.assertEqual(len(listed), 64)
        self.assertEqual(len(set(listed)), 64)

    def test_the_amino_acids_agree_with_biopython(self):
        """The table says which codon makes which amino acid. So does Biopython."""
        for amino, order in codons.ranked().items():
            for codon in order:
                self.assertEqual(silent.protein(codon), amino,
                                 f"{codon} is filed under {amino}")

    def test_the_yeast_favourites_are_right(self):
        """The test that fails if the package hands back another organism.

        Genome-wide yeast usage, which is not the same as the preference of
        highly expressed yeast genes: those favour AAG for lysine, the whole
        genome favours AAA. This asserts what the table actually says.
        """
        best = {amino: order[0] for amino, order in codons.ranked().items()}
        self.assertEqual(best["L"], "TTG")
        self.assertEqual(best["R"], "AGA")
        self.assertEqual(best["E"], "GAA")
        self.assertEqual(best["K"], "AAA")

    def test_ties_break_alphabetically(self):
        order = codons.ranked(9606)["R"]
        self.assertEqual(order[:2], ["AGA", "AGG"])  # both 0.21 in the table

    def test_the_method_string_names_the_table(self):
        self.assertIn(str(codons.YEAST), codons.method())
        self.assertIn("python_codon_tables-", codons.method())


class TestTheRewrite(unittest.TestCase):
    def setUp(self):
        self.genes = sg.read_genes(GENES)

    def test_the_protein_is_identical(self):
        for name, sequence, _ in self.genes:
            new, _ = codons.recode(sequence)
            self.assertEqual(str(Seq(new).translate()),
                             str(Seq(sequence).translate()), name)

    def test_the_length_is_identical(self):
        for name, sequence, _ in self.genes:
            new, _ = codons.recode(sequence)
            self.assertEqual(len(new), len(sequence), name)

    def test_the_gene_keeps_its_own_stop_codon(self):
        """Not merely some stop codon. Biopython's optimize() swaps it; we do not."""
        for name, sequence, _ in self.genes:
            new, _ = codons.recode(sequence)
            self.assertEqual(new[-3:], sequence[-3:], name)
            self.assertEqual(new[:3], "ATG", name)

    def test_the_result_is_still_a_type_3_cds(self):
        for name, sequence, _ in self.genes:
            new, _ = codons.recode(sequence)
            self.assertEqual(cds.problems(new, "full"), [], name)

    def test_the_cai_goes_up(self):
        for name, sequence, _ in self.genes:
            new, _ = codons.recode(sequence)
            self.assertGreater(metrics.cai(new), metrics.cai(sequence), name)

    def test_it_agrees_with_biopython_on_the_interior(self):
        """Pins the borrowed semantics.

        Only on a peptide with no amino acid three times in a row, because the
        run cap in lib/codons.py deliberately differs from plain argmax there,
        and only on the interior, because Biopython replaces the stop codon.
        """
        peptide = "ACDEFGHIKLMNPQRSTVWY" * 3
        ours, _ = codons.recode(gene(peptide))
        theirs = str(full_index().optimize(gene(peptide)))
        self.assertEqual(ours[3:-3], theirs[3:-3])

    def test_the_same_gene_always_gives_the_same_bases(self):
        for _, sequence, _ in self.genes:
            self.assertEqual(codons.recode(sequence)[0],
                             codons.recode(sequence)[0])

    def test_the_golden_sequence_has_not_moved(self):
        self.assertEqual(codons.recode(GOLDEN_IN)[0], GOLDEN_OUT)

    def test_the_table_really_drives_the_result(self):
        yeast, _ = codons.recode(GOLDEN_IN, codons.YEAST)
        human, _ = codons.recode(GOLDEN_IN, 9606)
        self.assertNotEqual(yeast, human)

    def test_a_run_of_one_amino_acid_is_not_one_repeated_codon(self):
        new, _ = codons.recode("ATG" + "CAG" * 10 + "TAA")
        self.assertNotIn("CAACAACAA", new)
        self.assertEqual(str(Seq(new).translate()), "M" + "Q" * 10 + "*")

    def test_the_rewrite_can_create_a_cut_site(self):
        """Why ytk-remove-cut-sites has to run after this, and is not optional.

        Trp-Ser-Gln recodes to TGGTCTCAA, which holds GGTCTC. Six three-codon
        windows of yeast top codons do this, all of them BsaI. A real case, not a
        theoretical one, so the step after this is mandatory.
        """
        clean = "ATG" + "TGGAGCCAG" + "GCTGCTGCT" + "TAA"
        self.assertEqual([], enzymes.in_sequence(clean))
        new, _ = codons.recode(clean)
        self.assertEqual(["BsaI"], enzymes.in_sequence(new))
        # And the step after it clears that again, without touching the protein.
        fixed, _, left = silent.remove_sites(new)
        self.assertEqual([], left)
        self.assertEqual(str(Seq(fixed).translate()), str(Seq(clean).translate()))

    def test_changes_are_reported_per_codon(self):
        _, changes = codons.recode(GOLDEN_IN)
        self.assertTrue(changes)
        for change in changes:
            self.assertNotEqual(change["was"], change["now"])
            self.assertEqual(silent.protein(change["was"]), change["amino_acid"])
            self.assertEqual(silent.protein(change["now"]), change["amino_acid"])


class TestWhatItRefuses(unittest.TestCase):
    def test_a_frameshifted_sequence(self):
        with self.assertRaises(ValueError):
            codons.recode("ATGAAAT")

    def test_letters_that_are_not_acgt(self):
        with self.assertRaises(ValueError):
            codons.recode("ATGNNNTAA")

    def test_no_start_codon(self):
        with self.assertRaises(ValueError):
            codons.recode("GGGAAATAA")

    def test_a_stop_codon_in_the_middle(self):
        """A warning in ytk-cds-qc, a refusal here.

        A sequence in the wrong frame can still be a whole number of codons that
        starts ATG and ends in a stop, and it rewrites into something that looks
        clean and means nothing.
        """
        with self.assertRaises(ValueError):
            codons.recode("ATGAAATAAAAATAA")


class TestTheMeasurements(unittest.TestCase):
    def test_a_planted_repeat_is_found_with_both_positions(self):
        unit = "AAGCTTGCATGCCTGCAGGT"  # 20 bases
        length, first, second = metrics.longest_repeat("GG" + unit + unit + "CC")
        self.assertEqual(length, 20)
        self.assertEqual((first, second), (3, 23))

    def test_nothing_short_of_the_threshold_counts_as_a_repeat(self):
        unit = "ACGTACGTAC"  # 10 bases, under REPEAT
        self.assertEqual(metrics.longest_repeat(unit + "TTTT" + unit)[0], 0)

    def test_a_planted_run_is_found(self):
        self.assertEqual(metrics.longest_run("GCGC" + "A" * 9 + "GCGC"),
                         ("A", 9, 5))

    def test_gc_windows_find_a_low_stretch_a_whole_gene_average_hides(self):
        sequence = "GC" * 60 + "AT" * 60
        low, high = metrics.gc_windows(sequence)
        self.assertEqual((low, high), (0.0, 100.0))
        self.assertEqual(metrics.gc(sequence), 50.0)

    def test_top_codons_score_one(self):
        """Also checks the index is built right, per organism."""
        for taxid in (codons.YEAST, 316407, 9606):
            best = codons._choices(taxid)
            sequence = "ATG" + "".join(best[a][0] for a in "ACDEFGHIKLNPQRSTVY") + "TAA"
            self.assertEqual(metrics.cai(sequence, taxid), 1.0, taxid)

    def test_the_stop_codon_does_not_change_the_score(self):
        """Sharp & Li exclude it, and a gene cannot choose a different one."""
        scores = {metrics.cai("ATGCTGGAGAAG" + stop)
                  for stop in ("TAA", "TAG", "TGA")}
        self.assertEqual(len(scores), 1)

    def test_no_cai_for_a_sequence_it_cannot_score(self):
        self.assertIsNone(metrics.cai("ATGAA"))        # not whole codons
        self.assertIsNone(metrics.cai("ATGNNNTAA"))    # not ACGT
        self.assertIsNone(metrics.cai("ATGTAA"))       # nothing scorable in it
        self.assertIsNone(metrics.cai(""))

    def test_measure_fills_every_column(self):
        measured = metrics.measure(GOLDEN_IN)
        self.assertEqual(sorted(measured), sorted(metrics.COLUMNS))


if __name__ == "__main__":
    unittest.main(verbosity=2)
