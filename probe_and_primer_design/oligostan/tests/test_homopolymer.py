import unittest
from unittest.mock import patch

import pandas as pd
from Bio.Seq import Seq

from oligostan.audit_run import audit_homopolymers
from oligostan.filters import is_ok_4_homopolymer, longest_homopolymer_run
from oligostan.mouse_smifish import generate_candidates


class HomopolymerTests(unittest.TestCase):
    def test_four_passes_and_five_fails_for_each_base(self):
        for base in "ATCG":
            with self.subTest(base=base):
                flank = next(other for other in "ATCG" if other != base)
                self.assertTrue(is_ok_4_homopolymer(flank + base * 4 + flank))
                self.assertFalse(is_ok_4_homopolymer(flank + base.lower() * 5 + flank))
                self.assertEqual(longest_homopolymer_run(flank + base.lower() * 5), 5)

    def test_candidate_generation_excludes_homopolymers_before_blast(self):
        target = "ACGT" * 13 + "A" * 26
        antisense = str(Seq(target).reverse_complement())
        probes = [(26, 1.0, pos, antisense[pos - 1:pos + 25]) for pos in (1, 30)]
        metric = {"GCpc": 0.5, "GCFilter": 1, "PNASFilter": 1,
                  "aCompFilter": 1, "aStackFilter": 1, "cCompFilter": 1,
                  "cStackFilter": 1, "cSpecStackFilter": 1, "NbOfPNAS": 5,
                  "dG37": -32.0}
        relaxed_metric = {**metric, "PNASFilter": 0, "aCompFilter": 0,
                          "aStackFilter": 0, "cStackFilter": 0, "NbOfPNAS": 2}
        model = {"Spen": {"gene": "Spen", "gene_id": "Spen", "transcript_id": "t1",
                          "chrom": "chr1", "strand": "+", "exon": [(1, len(target))],
                          "intron": []}}
        with patch("oligostan.mouse_smifish.FastaIndex") as fasta, \
             patch("oligostan.mouse_smifish.tile_intervals", return_value=[(1, len(target))]), \
             patch("oligostan.mouse_smifish.get_probes_from_rna_dg37", return_value=probes), \
             patch("oligostan.mouse_smifish.process_probes_for_output",
                   return_value=[relaxed_metric, metric]):
            fasta.return_value.fetch.return_value = target
            candidates = generate_candidates(model, fasta="unused")
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates.probe_seq.iloc[0], probes[1][3])
        self.assertLessEqual(candidates.max_homopolymer_run.iloc[0], 4)

    def test_audit_rejects_failing_or_misreported_run(self):
        good = pd.DataFrame([{"probe_id": "q1", "probe_seq": "ATCGATCG", "max_homopolymer_run": 1}])
        audit_homopolymers(good, 4, "fixture")
        bad = pd.DataFrame([{"probe_id": "q2", "probe_seq": "ATCGAAAAA", "max_homopolymer_run": 5}])
        with self.assertRaisesRegex(AssertionError, "q2"):
            audit_homopolymers(bad, 4, "fixture")
        good.loc[0, "max_homopolymer_run"] = 2
        with self.assertRaisesRegex(AssertionError, "q1"):
            audit_homopolymers(good, 4, "fixture")


if __name__ == "__main__":
    unittest.main()
