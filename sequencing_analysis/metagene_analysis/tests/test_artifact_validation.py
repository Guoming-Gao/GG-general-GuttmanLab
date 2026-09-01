import json
from pathlib import Path

import numpy as np
import pandas as pd

from spen_metagene.artifact_validation import (
    assign_tertiles,
    feature_intervals,
    feature_lengths,
    filter_ca,
    peak_count_universe,
    peak_statistics,
)


def _model(strand="+"):
    return {
        "transcript_id": "tx1", "gene_id": "g1", "gene_name": "Gene1",
        "chrom": "chr1", "strand": strand, "start": 0, "end": 300,
        "gene_type": "protein_coding", "transcript_type": "protein_coding",
        "exons": [[0, 100], [200, 300]], "cds": [[20, 80], [220, 280]],
    }


def test_feature_geometry_respects_transcript_strand():
    plus = feature_intervals(_model("+")); minus = feature_intervals(_model("-"))
    assert plus["first_exon"] == [[0, 100]]
    assert minus["first_exon"] == [[200, 300]]
    assert plus["five_prime_utr"] == [[0, 20]]
    assert minus["five_prime_utr"] == [[280, 300]]
    assert plus["all_introns"] == minus["all_introns"] == [[100, 200]]


def test_ca_filters_are_strict_and_deduplicate_windows(tmp_path):
    columns = ["window", "strand", "feature name", "sample count", "input count",
               "enrichment (intra)", "enrichment (inter)", "p-val (intra)", "p-val (inter)"]
    rows = [
        ["chr1:0-100", "+", "x", 5, 1, 2.01, 2.01, .009, .009],
        ["chr1:0-100", "+", "x", 6, 1, 3, 3, .001, .001],
        ["chr1:100-200", "+", "x", 5, 1, 2.0, 3, .001, .001],
        ["chr1:200-300", "+", "x", 5, 1, 3, 3, .01, .001],
    ]
    p = tmp_path / "ca.tsv"; pd.DataFrame(rows, columns=columns).to_csv(p, sep="\t", index=False)
    out = filter_ca(p)
    assert out.window.tolist() == ["chr1:0-100"]


def test_absolute_peak_count_keeps_zeros_and_is_not_length_normalized():
    models = {"g1": _model(), "g2": {**_model(), "gene_id": "g2", "gene_name": "Gene2", "start": 400, "end": 900,
                                            "exons": [[400, 500], [800, 900]], "cds": [[420, 480], [820, 880]]}}
    lengths = feature_lengths(models)
    expr = pd.DataFrame({"gene_id": ["g1", "g2"], "expression": [1, 10],
                         "log10_expression_plus_1": np.log10([2, 11]),
                         "expression_tertile": ["Low", "High"], "expression_source": ["x", "x"]})
    assigned = pd.DataFrame({"gene_id": ["g1", "g1"], "feature": ["first_exon", "first_exon"],
                             "window": ["chr1:0-100", "chr1:100-200"], "callset": ["pooled", "pooled"]})
    out = peak_count_universe(assigned, lengths, expr, "pooled")
    assert out[(out.gene_id == "g1") & (out.feature == "first_exon")].absolute_CA_peak_count.iloc[0] == 2
    assert out[(out.gene_id == "g2") & (out.feature == "first_exon")].absolute_CA_peak_count.iloc[0] == 0


def test_tertile_ties_are_not_split():
    d = pd.DataFrame({"gene_id": list("abcdef"), "expression": [1, 1, 2, 2, 3, 3]})
    out, _ = assign_tertiles(d)
    assert out.groupby("expression").expression_tertile.nunique().max() == 1


def test_ca_association_statistics_exclude_zero_peak_features_but_report_them_as_qc():
    n_positive = 12
    d = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(n_positive + 20)],
        "feature": "all_introns",
        "feature_length": 1000,
        "log10_expression_plus_1": np.arange(n_positive + 20, dtype=float),
        "expression_tertile": (["Low"] * 4 + ["Middle"] * 4 + ["High"] * 4 + ["Low"] * 20),
        "absolute_CA_peak_count": list(range(1, n_positive + 1)) + [0] * 20,
        "any_CA_peak": [True] * n_positive + [False] * 20,
    })
    summary, models, _ = peak_statistics(d, bootstrap=10, seed=7)
    row = summary[(summary.feature == "all_introns") & (summary.tertile == "Low")].iloc[0]
    assert row.n_peak_positive_genes == 4
    assert row.n_all_eligible_genes == 24
    assert row.fraction_zero_QC == 20 / 24
    assert summary[summary.feature == "all_introns"].spearman_rho_peak_positive.iloc[0] == 1
    assert set(models.endpoint.dropna()) == {"positive_peak_count_negative_binomial"}
