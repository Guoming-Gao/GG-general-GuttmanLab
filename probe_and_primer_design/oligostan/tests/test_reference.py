import gzip
import tempfile
import unittest
from pathlib import Path

from oligostan.mouse_smifish import FastaIndex, load_models, quality_tier
from oligostan.filters import is_it_ok_4_c_spec_stack


class ReferenceInputTests(unittest.TestCase):
    def test_quality_tiers_relax_project_rules_in_order(self):
        metric = {"GCpc": 0.5, "aCompFilter": 1, "aStackFilter": 1,
                  "cCompFilter": 0, "cStackFilter": 1, "cSpecStackFilter": 0}
        self.assertEqual(quality_tier(metric), 0)
        metric["aCompFilter"] = 0
        self.assertEqual(quality_tier(metric), 1)
        metric["aStackFilter"] = 0
        self.assertEqual(quality_tier(metric), 2)
        metric["cStackFilter"] = 0
        self.assertEqual(quality_tier(metric), 3)

    def test_r_fifth_rule_only_checks_first_six_windows(self):
        self.assertTrue(is_it_ok_4_c_spec_stack("ATATATATATATCCCCATATATATATAT"))

    def test_fasta_index_crosses_wrapped_lines_without_coordinate_shift(self):
        with tempfile.TemporaryDirectory() as tmp:
            fasta = Path(tmp) / "tiny.fa"
            fasta.write_bytes(b">chr1\nACGT\nTGCA\nAT\n")
            Path(str(fasta) + ".fai").write_text("chr1\t10\t6\t4\t5\n")
            indexed = FastaIndex(fasta)
            self.assertEqual(indexed.fetch("chr1", 1, 10), "ACGTTGCAAT")
            self.assertEqual(indexed.fetch("chr1", 4, 6), "TTG")
            self.assertEqual(indexed.fetch("chr1", 9, 10), "AT")

    def test_alternative_exon_is_excluded_from_intron_set(self):
        with tempfile.TemporaryDirectory() as tmp:
            gtf = Path(tmp) / "tiny.gtf.gz"
            lines = []
            for tid, spans in (("t1", ((1, 40), (201, 250))),
                               ("t2", ((1, 30), (90, 120), (201, 220)))):
                for start, end in spans:
                    lines.append(f'chr1\trefGene\texon\t{start}\t{end}\t.\t+\t.\tgene_id "G"; transcript_id "{tid}"; gene_name "Spen";\n')
            with gzip.open(gtf, "wt") as out:
                out.writelines(lines)
            model = load_models(gtf, genes=("Spen",))["Spen"]
            self.assertEqual(model["transcript_id"], "t1")
            self.assertEqual(model["intron"], [(41, 89), (121, 200)])


if __name__ == "__main__":
    unittest.main()
