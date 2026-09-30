"""BLAST-first minimum-span selection, an add-on to the Oligostan port.

The input is the complete BLAST-tested candidate table. Oligostan's tiled
design produces nonoverlapping probes within each gene/region, so every group
of 30 consecutive passing probes is a valid window to compare.
"""

import pandas as pd


def expected_sets(models):
    """Return annotated gene/region pairs, omitting absent introns."""
    return [(gene, region) for gene, model in models.items()
            for region in ("exon", "intron") if model[region]]


def select_minimum_span_sets(verified, models, probes_per_set=30):
    """Select the narrowest genomic window from *all* BLAST-passing probes.

    No preselection by score, quality tier, or old 40-probe limit is applied.
    Quality tiers and scores break ties only after genomic span.
    """
    if probes_per_set < 1:
        raise ValueError("probes_per_set must be positive")
    if verified.probe_id.duplicated().any():
        raise ValueError("Probe IDs must be unique")
    selected_parts, summary_rows = [], []
    for gene, region in expected_sets(models):
        tested = verified[(verified.gene == gene) & (verified.region == region)]
        pool = tested[tested.blast_verified].sort_values(
            ["chrom", "start", "end", "probe_id"], kind="stable").reset_index(drop=True)
        if pool.chrom.nunique() > 1:
            raise ValueError(f"Multiple chromosomes for {gene} {region}")
        if len(pool) > 1 and (pool.start.iloc[1:].to_numpy() <= pool.end.iloc[:-1].to_numpy()).any():
            raise ValueError(f"Overlapping Oligostan candidates for {gene} {region}")
        if len(pool) >= probes_per_set:
            windows = []
            for i in range(len(pool) - probes_per_set + 1):
                subset = pool.iloc[i:i + probes_per_set]
                windows.append((
                    int(subset.end.iloc[-1] - subset.start.iloc[0] + 1),
                    int(subset.quality_tier.max()), int(subset.quality_tier.sum()),
                    -float(subset.dGScore.sum()), int(subset.start.iloc[0]), i,
                ))
            chosen = pool.iloc[min(windows)[-1]:min(windows)[-1] + probes_per_set].copy()
        else:
            chosen = pool.copy()
        selected_parts.append(chosen)
        summary_rows.append({
            "gene": gene, "region": region,
            "blast_tested_count": len(tested), "blast_verified_count": len(pool),
            "selected_count": len(chosen), "minimum_requested": probes_per_set,
            "meets_minimum": len(chosen) == probes_per_set,
            "chrom": chosen.chrom.iloc[0] if len(chosen) else models[gene]["chrom"],
            "span_start": int(chosen.start.min()) if len(chosen) else None,
            "span_end": int(chosen.end.max()) if len(chosen) else None,
            "genomic_span_bp": int(chosen.end.max() - chosen.start.min() + 1) if len(chosen) else None,
            "max_quality_tier": int(chosen.quality_tier.max()) if len(chosen) else None,
        })
    selected = pd.concat(selected_parts, ignore_index=True) if selected_parts else verified.iloc[:0].copy()
    selected["set_id"] = selected.gene + "_" + selected.region + "_mm10"
    return selected, pd.DataFrame(summary_rows)
