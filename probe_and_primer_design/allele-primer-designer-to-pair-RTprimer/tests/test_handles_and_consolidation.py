from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from design_allele_primers_to_pair_rtprimer import DEFAULT_CONFIG
from utils.designer import HANDLE_SEQUENCES, _expand_handle
from utils.report import write_consolidated_html_report


class HandleSelectionTests(unittest.TestCase):
    def setUp(self) -> None:
        self.insert = {"Primer_Insert_Seq": "ACGTACGT"}

    def test_default_handle_is_2puni(self) -> None:
        self.assertEqual(DEFAULT_CONFIG["handle"], "2PUNI")
        rows = _expand_handle([self.insert], DEFAULT_CONFIG["handle"])
        self.assertEqual([row["Handle_Name"] for row in rows], ["2PUNI"])

    def test_2pbc_handle_is_selected_exclusively(self) -> None:
        rows = _expand_handle([self.insert], "2PBC")
        self.assertEqual([row["Handle_Name"] for row in rows], ["2PBC"])
        self.assertEqual(rows[0]["Final_Order_Sequence"], HANDLE_SEQUENCES["2PBC"] + "ACGTACGT")

    def test_unknown_handle_is_rejected(self) -> None:
        with self.assertRaises(ValueError):
            _expand_handle([self.insert], "unknown")


class ConsolidatedReportTests(unittest.TestCase):
    def test_report_contains_successful_oligo_and_failed_gene(self) -> None:
        batch_rows = [
            {"Gene": "Xist", "Status": "success", "Selected_Transcript": "NR_1", "Strand": "-", "Informative_SNPs": 2, "Candidate_Inserts": 3, "Top_Inserts": 1, "Top_Oligos": 1, "Error": ""},
            {"Gene": "Tsix", "Status": "failed", "Selected_Transcript": "", "Strand": "", "Informative_SNPs": 0, "Candidate_Inserts": 0, "Top_Inserts": 0, "Top_Oligos": 0, "Error": "test failure"},
        ]
        oligos = [{"Gene": "Xist", "Transcript_ID": "NR_1", "Strand": "-", "Insert_Rank": 1, "Primer_Insert_Seq": "ACGT", "Handle_Name": "2PUNI", "Final_Order_Sequence": HANDLE_SEQUENCES["2PUNI"] + "ACGT", "SNP_cDNA": 10, "SNP_Genomic": 20, "SNP_Distance_To_RNA_3p": 30, "Primer_To_SNP_Distance_Toward_PolyDT": 5, "Primer_Tm": 60, "Primer_GC": 50, "Primer_Specificity_Pass": True}]
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "report.html"
            write_consolidated_html_report(output, batch_rows, oligos, {"handle": "2PUNI", "snp_window_toward_polydt": 100, "blast_specificity": True})
            html = output.read_text()
        self.assertIn("Xist", html)
        self.assertIn("Tsix", html)
        self.assertIn("test failure", html)
        self.assertIn(HANDLE_SEQUENCES["2PUNI"] + "ACGT", html)


if __name__ == "__main__":
    unittest.main()
