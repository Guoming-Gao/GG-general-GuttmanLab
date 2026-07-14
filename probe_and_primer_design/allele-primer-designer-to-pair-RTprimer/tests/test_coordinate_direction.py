from __future__ import annotations

import sys
import unittest
from pathlib import Path


PACKAGE_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PACKAGE_DIR))

from utils.coordinate_rules import primer_to_snp_distance_toward_polydt, snp_distance_to_rna_3p
from utils.models import CdnaExonSegment, Exon, TranscriptTemplate
from utils.transcripts import cdna_interval_to_genomic_segments, map_cdna_pos_to_genomic, map_genomic_pos_to_cdna


class CoordinateDirectionTests(unittest.TestCase):
    def test_plus_strand_cdna_increases_with_genome(self) -> None:
        exon = Exon("chr1", 100, 199, 1)
        template = TranscriptTemplate("GeneP", "txP", "chr1", "+", "A" * 100, (CdnaExonSegment(exon, 1, 100),))

        self.assertEqual(map_cdna_pos_to_genomic(template, 1), 100)
        self.assertEqual(map_cdna_pos_to_genomic(template, 100), 199)
        self.assertEqual(map_genomic_pos_to_cdna(template, 150), 51)
        self.assertEqual(cdna_interval_to_genomic_segments(template, 20, 40), [("chr1", 119, 139, "+")])

    def test_minus_strand_cdna_increases_while_genome_decreases(self) -> None:
        exon = Exon("chrX", 100, 199, 1)
        template = TranscriptTemplate("GeneM", "txM", "chrX", "-", "A" * 100, (CdnaExonSegment(exon, 1, 100),))

        self.assertEqual(map_cdna_pos_to_genomic(template, 1), 199)
        self.assertEqual(map_cdna_pos_to_genomic(template, 100), 100)
        self.assertEqual(map_genomic_pos_to_cdna(template, 150), 50)
        self.assertEqual(cdna_interval_to_genomic_segments(template, 20, 40), [("chrX", 160, 180, "-")])

    def test_snp_must_be_toward_polydt_from_primer_3p_end(self) -> None:
        self.assertEqual(primer_to_snp_distance_toward_polydt(40, 50, 100), 10)
        self.assertEqual(primer_to_snp_distance_toward_polydt(40, 140, 100), 100)
        self.assertIsNone(primer_to_snp_distance_toward_polydt(40, 141, 100))
        self.assertIsNone(primer_to_snp_distance_toward_polydt(60, 50, 100))
        self.assertIsNone(primer_to_snp_distance_toward_polydt(50, 50, 100))

    def test_snp_distance_to_rna_3p_is_spliced_distance(self) -> None:
        self.assertEqual(snp_distance_to_rna_3p(1000, 1000), 0)
        self.assertEqual(snp_distance_to_rna_3p(1000, 900), 100)


if __name__ == "__main__":
    unittest.main()
