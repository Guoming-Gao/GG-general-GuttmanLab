from pathlib import Path

from spen_metagene.annotation import parse_representative_transcripts
from spen_metagene.clap import feature_membership


def test_coding_preference_and_strand_features(tmp_path: Path):
    g=tmp_path/"x.gtf"
    g.write_text(
      'chr1\tx\ttranscript\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "non"; gene_name "G1"; transcript_type "lncRNA";\n'
      'chr1\tx\texon\t101\t400\t.\t+\t.\tgene_id "g1"; transcript_id "non";\n'
      'chr1\tx\ttranscript\t101\t300\t.\t+\t.\tgene_id "g1"; transcript_id "coding"; gene_name "G1"; transcript_type "protein_coding";\n'
      'chr1\tx\texon\t101\t150\t.\t+\t.\tgene_id "g1"; transcript_id "coding";\n'
      'chr1\tx\texon\t201\t300\t.\t+\t.\tgene_id "g1"; transcript_id "coding";\n'
      'chr1\tx\tCDS\t121\t250\t.\t+\t0\tgene_id "g1"; transcript_id "coding";\n')
    t=parse_representative_transcripts(str(g))["g1"]
    assert t.transcript_id=="coding"
    td=t.__dict__
    assert "five_prime_utr" in feature_membership(td,110)
    assert "first_intron" in feature_membership(td,175)
    assert "three_prime_utr" in feature_membership(td,275)


def test_minus_strand_feature_order_is_transcription_oriented():
    t={"strand":"-","exons":[[100,150],[200,250],[300,350]],"cds":[[120,330]]}
    assert "first_exon" in feature_membership(t,325)
    assert "first_intron" in feature_membership(t,275)
    assert "later_introns" in feature_membership(t,175)
    assert "five_prime_utr" in feature_membership(t,340)
    assert "three_prime_utr" in feature_membership(t,110)
