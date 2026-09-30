import unittest
from unittest.mock import patch

import pandas as pd

from oligostan.minimum_span import expected_sets, select_minimum_span_sets
from oligostan.mouse_smifish import adaptive_blast_verify


class MinimumSpanTests(unittest.TestCase):
    def setUp(self):
        self.models = {"Spen": {"exon": [(1, 20000)], "intron": [], "chrom": "chr1"},
                       "Malat1": {"exon": [(1, 7000)], "intron": [], "chrom": "chr19"}}

    def frame(self, starts, passing=True):
        return pd.DataFrame([{
            "probe_id": f"q{i:04d}", "gene": "Spen", "region": "exon",
            "chrom": "chr1", "start": start, "end": start + 19,
            "blast_verified": passing, "quality_tier": 0, "dGScore": 1.0,
        } for i, start in enumerate(starts)])

    def test_full_passing_pool_can_choose_probes_beyond_old_forty_cap(self):
        wide = list(range(100, 3100, 100))
        compact = list(range(10000, 10900, 30))
        base = self.frame(wide)
        _, before = select_minimum_span_sets(base, {"Spen": self.models["Spen"]})
        full = self.frame(wide + compact)
        chosen, after = select_minimum_span_sets(full, {"Spen": self.models["Spen"]})
        self.assertEqual(len(chosen), 30)
        self.assertTrue(all(chosen.start >= 10000))
        self.assertLess(after.genomic_span_bp.iloc[0], before.genomic_span_bp.iloc[0])
        extended = self.frame(wide + compact + [15000])
        _, unchanged = select_minimum_span_sets(extended, {"Spen": self.models["Spen"]})
        self.assertEqual(after.genomic_span_bp.iloc[0], unchanged.genomic_span_bp.iloc[0])

    def test_blast_failing_candidates_cannot_enter_window(self):
        passing = self.frame(list(range(100, 3100, 100)))
        failed = self.frame(list(range(10000, 10900, 30)), passing=False)
        failed.probe_id = [f"f{i:04d}" for i in range(len(failed))]
        chosen, _ = select_minimum_span_sets(pd.concat([passing, failed]),
                                              {"Spen": self.models["Spen"]})
        self.assertTrue(chosen.blast_verified.all())
        self.assertTrue(all(chosen.start < 10000))

    def test_no_malat1_intron_set(self):
        self.assertEqual(expected_sets(self.models), [("Spen", "exon"), ("Malat1", "exon")])

    def test_equal_span_ties_use_tier_then_score(self):
        rows = self.frame([100, 130, 160])
        rows.loc[0, "quality_tier"] = 1
        chosen, _ = select_minimum_span_sets(rows, {"Spen": self.models["Spen"]}, probes_per_set=2)
        self.assertEqual(chosen.start.tolist(), [130, 160])
        rows.quality_tier = 0
        rows.dGScore = [1.0, 1.0, 3.0]
        chosen, _ = select_minimum_span_sets(rows, {"Spen": self.models["Spen"]}, probes_per_set=2)
        self.assertEqual(chosen.start.tolist(), [130, 160])

    def test_batched_blast_visits_every_candidate(self):
        rows = self.frame(list(range(100, 36100, 30)))
        with patch("oligostan.mouse_smifish.blast_verify") as blast:
            blast.side_effect = lambda batch, **kwargs: (batch.copy(), pd.DataFrame(), ["blastn"])
            verified, _, commands = adaptive_blast_verify(rows, batch_size=500)
        self.assertEqual(len(verified), len(rows))
        self.assertEqual(len(commands), 3)
        self.assertEqual(blast.call_count, 3)


if __name__ == "__main__":
    unittest.main()
