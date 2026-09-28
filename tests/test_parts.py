"""Tests for the data tables and the entry vector.

These check the repo against the published toolkit rather than against itself:
the 96 plasmid maps in reference/ytk_plasmids/ are the ground truth.

    python3 tests/test_parts.py
"""
import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "lib"))
import enzymes    # noqa: E402
import flanks     # noqa: E402
import snapgene as sg  # noqa: E402

PLASMIDS = ROOT / "reference" / "ytk_plasmids"
ENTRY_VECTOR = PLASMIDS / "pYTK001.gb"

# Measured from pYTK001.gb: 2676 bp total, a 1030 bp GFP dropout that the part
# replaces, and 1646 bp kept. Those are the contributions to the finished
# circle. pydna counts a cut piece's sticky ends as part of it, so a cut piece
# measures four bases more than it contributes, and two pieces joined lose eight.
ENTRY_VECTOR_LENGTH = 2676
BACKBONE_LENGTH = 1646
DROPOUT_LENGTH = 1030
OVERHANG = 4


def table(name):
    with open(ROOT / "data" / name, newline="") as fh:
        return [row for row in csv.DictReader(fh, delimiter="\t")
                if not row["part_type"].startswith("#")]


class TestTheOverhangTable(unittest.TestCase):
    """data/ytk_overhangs.tsv against data/ytk_part_types.tsv, which came
    straight from the paper."""

    def setUp(self):
        self.adapters = {row["part_type"]: row for row in table("ytk_overhangs.tsv")}
        self.published = {row["part_type"]: row for row in table("ytk_part_types.tsv")}

    def test_every_adapter_carries_the_published_overhang(self):
        """An adapter may add extra bases the paper asks for, such as the BglII
        site on Type 2 or the NotI sites on Type 8, but the four overhang bases
        have to be exactly right or the part goes in the wrong slot.

        The overhang sits at the junction, so it comes first in a left adapter
        and last in a right one, with any extra bases on the inward side.
        """
        for part, published in self.published.items():
            if part not in self.adapters:
                continue
            with self.subTest(part=part):
                left = self.adapters[part]["left_adapter"].upper()
                right = self.adapters[part]["right_adapter"].upper()
                upstream = published["upstream"].upper()
                downstream = published["downstream"].upper()
                # Types 3 and 3a store a bare T, because the sequence's own ATG
                # finishes the TATG overhang.
                if left != "T":
                    self.assertTrue(left.startswith(upstream),
                                    f"{left} should start with {upstream}")
                else:
                    self.assertEqual(upstream, "TATG")
                self.assertTrue(right.endswith(downstream),
                                f"{right} should end with {downstream}")

    def test_the_overhangs_chain_into_a_circle(self):
        chain = ["1", "2", "3", "4", "5", "6", "7", "8"]
        for part, following in zip(chain, chain[1:] + chain[:1]):
            with self.subTest(part=part):
                self.assertEqual(self.published[part]["downstream"],
                                 self.published[following]["upstream"])

    def test_both_tables_use_the_same_part_type_names(self):
        parts = {row["part_type"] for row in table("ytk_parts.tsv")}
        adapters = set(self.adapters) - {"custom"}
        self.assertTrue(adapters <= parts | {"3"},
                        f"only in the adapter table: {sorted(adapters - parts)}")

    def test_nothing_is_flagged_unsure_any_more(self):
        """Types 8 and 8a were marked unverified. The paper and pYTK083-085 and
        pYTK089-091 confirm both.

        No code reads this column any more - see NOTES.md. This test is what
        keeps the verification on record, so refilling the cell fails here
        rather than being silently ignored at run time.
        """
        self.assertEqual([row["part_type"] for row in table("ytk_overhangs.tsv")
                          if row.get("unsure")], [])


class TestNoSecondStopCodon(unittest.TestCase):
    """Plasmid_Generator.xlsx stores the Type 3 and 3b right pieces as TAGATCC.
    The TAG is an extra stop codon added whether or not the CDS has one."""

    def setUp(self):
        self.adapters = {row["part_type"]: row for row in table("ytk_overhangs.tsv")}

    def test_type_3_and_3b_right_adapters_are_plain_atcc(self):
        for part in ("3", "3b"):
            with self.subTest(part=part):
                self.assertEqual(self.adapters[part]["right_adapter"], "ATCC")

    def test_no_adapter_anywhere_holds_a_tag_before_atcc(self):
        for part, row in self.adapters.items():
            with self.subTest(part=part):
                self.assertNotIn("TAGATCC", row["right_adapter"].upper())

    def test_a_flanked_cds_keeps_exactly_its_own_stop_codon(self):
        """Searching the fragment for TAGATCC is not the test: a CDS whose own
        stop codon is TAG produces that string legitimately. What matters is
        that the insert between the adapters is the sequence as given, with
        nothing added.
        """
        for stop in ("TAA", "TAG", "TGA"):
            with self.subTest(stop=stop):
                gene = "ATGAAACGT" + stop
                fragment = flanks.flank(gene, self.adapters["3"]).upper()
                start = flanks.insert_offset(self.adapters["3"])
                self.assertEqual(fragment[start:start + len(gene)], gene)
                self.assertEqual(fragment[start + len(gene):start + len(gene) + 4],
                                 "ATCC")

    def test_the_spreadsheets_extra_tag_would_be_caught(self):
        """What the bug looked like: a TAG between the stop codon and ATCC, so
        two stop codons in a row."""
        broken = dict(self.adapters["3"], right_adapter="TAGATCC")
        fragment = flanks.flank("ATGAAACGTTAA", broken).upper()
        start = flanks.insert_offset(broken)
        after = fragment[start + len("ATGAAACGTTAA"):]
        self.assertTrue(after.startswith("TAG"))
        self.assertFalse(after.startswith("ATCC"))


class TestTheEntryVector(unittest.TestCase):
    """Numbers measured from pYTK001.gb, so a change to the vector or the
    cutting code shows up here."""

    @classmethod
    def setUpClass(cls):
        sys.path.insert(0, str(ROOT / "skills" / "ytk-clone" / "scripts"))
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "clone", ROOT / "skills" / "ytk-clone" / "scripts" / "clone.py")
        cls.clone = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.clone)
        cls.backbone = cls.clone.backbone_sequence()

    def test_the_entry_vector_is_2676_bp(self):
        self.assertEqual(len(self.backbone), ENTRY_VECTOR_LENGTH)

    def test_it_holds_two_bsmbi_sites_and_no_bsai_site(self):
        self.assertEqual(enzymes.count(self.backbone, "BsmBI"), 2)
        self.assertEqual(enzymes.count(self.backbone, "BsaI"), 0)

    def test_cutting_it_leaves_the_1646_bp_piece(self):
        kept = self.clone.cut_backbone(self.backbone)
        self.assertEqual(len(kept) - OVERHANG, BACKBONE_LENGTH)

    def test_the_two_pieces_account_for_the_whole_vector(self):
        pieces = self.clone.Dseqrecord(self.backbone, circular=True).cut(
            self.clone.BsmBI)
        self.assertEqual(sum(len(piece) for piece in pieces) - 2 * OVERHANG,
                         ENTRY_VECTOR_LENGTH)

    def test_the_dropout_is_the_rest(self):
        self.assertEqual(ENTRY_VECTOR_LENGTH - BACKBONE_LENGTH, DROPOUT_LENGTH)

    def test_a_part_plasmid_is_the_backbone_plus_the_insert(self):
        """A free check on every map: if the arithmetic does not hold, something
        was added or lost."""
        adapters = {row["part_type"]: row
                    for row in table("ytk_overhangs.tsv")}["3"]
        for gene in ("ATGAAACGTTAA", "ATG" + "AAACGT" * 20 + "TGA"):
            with self.subTest(length=len(gene)):
                fragment = flanks.flank(gene, adapters)
                insert = self.clone.cut_insert(fragment)
                plasmid = self.clone.assemble(fragment, self.backbone)
                self.assertEqual(len(plasmid),
                                 BACKBONE_LENGTH + len(insert) - OVERHANG)

    def test_pad_length_does_not_change_the_plasmid(self):
        """The pad sits outside both BsmBI sites, so it is cut off and thrown
        away. A fragment ordered with a 10-base pad and a PCR product with a
        4-base pad have to give the same plasmid."""
        adapters = {row["part_type"]: row
                    for row in table("ytk_overhangs.tsv")}["3"]
        gene = "ATG" + "AAACGTTTCGAC" * 8 + "TAA"
        built = []
        for pad in (flanks.MIN_PAD, 7, flanks.FULL_PAD):
            fragment = (flanks.forward_pad(pad) + flanks.FORWARD_SCAFFOLD
                        + adapters["left_adapter"] + gene
                        + adapters["right_adapter"] + flanks.RIGHT_SCAFFOLD
                        + flanks.reverse_complement(flanks.reverse_pad(pad)))
            built.append(self.clone.assemble(fragment, self.backbone))
        self.assertEqual(len(set(built)), 1, "pad length changed the plasmid")


class TestThePartsTableMatchesThePublishedPlasmids(unittest.TestCase):
    def test_all_ninety_six_plasmids_are_in_the_repo(self):
        self.assertEqual(len(list(PLASMIDS.glob("pYTK*.gb"))), 96)

    def test_the_table_rebuilds_byte_for_byte(self):
        """Rebuilds data/ytk_parts.tsv from the published maps and compares. If
        the generator and the table ever disagree, this is where it shows."""
        try:
            import pandas  # noqa: F401
        except ImportError:
            self.skipTest("needs requirements-dev.txt (pandas, xlrd)")
        committed = (ROOT / "data" / "ytk_parts.tsv").read_bytes()
        with tempfile.TemporaryDirectory() as tmp:
            copy = Path(tmp) / "ytk_parts.tsv"
            copy.write_bytes(committed)
            subprocess.run([sys.executable,
                            str(ROOT / "data" / "build_parts_table.py")],
                           check=True, capture_output=True)
            rebuilt = (ROOT / "data" / "ytk_parts.tsv").read_bytes()
            (ROOT / "data" / "ytk_parts.tsv").write_bytes(committed)
        self.assertEqual(rebuilt, committed)

    def test_the_table_has_a_row_for_every_plasmid(self):
        rows = table("ytk_parts.tsv")
        self.assertEqual(len(rows), 100)  # 96 parts + 4 entry-vector features

    def test_every_type_3_part_has_the_tatg_and_atcc_junctions(self):
        for row in table("ytk_parts.tsv"):
            if row["part_type"] == "3":
                with self.subTest(part=row["name"]):
                    self.assertEqual(row["junction_5"], "TATG")
                    self.assertEqual(row["junction_3"], "ATCC")


if __name__ == "__main__":
    unittest.main(verbosity=2)
