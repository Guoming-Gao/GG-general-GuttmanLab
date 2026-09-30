"""Genomic coverage figures for BLAST-verified Oligostan smiFISH sets."""

from pathlib import Path

from .mouse_smifish import progress


def write_coverage_report(root, models, selected, summary):
    """Write a summary CSV, one PNG per set, and a multipage PDF."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.ticker import FuncFormatter, MaxNLocator

    root = Path(root)
    plot_dir = root / "coverage_plots"
    plot_dir.mkdir(exist_ok=True)
    summary.to_csv(root / "coverage_summary.csv", index=False)
    pdf_path = root / "coverage_report.pdf"
    with PdfPages(pdf_path) as pdf, progress() as bar:
        task = bar.add_task("Drawing genomic coverage report", total=len(summary))
        for item in summary.itertuples(index=False):
            model = models[item.gene]
            probes = selected[(selected.gene == item.gene) & (selected.region == item.region)].sort_values("start")
            fig, (context, zoom) = plt.subplots(2, 1, figsize=(12, 6.4))
            color = "#2166ac" if item.region == "exon" else "#b35806"
            gene_start = min(start for start, _ in model["exon"])
            gene_end = max(end for _, end in model["exon"])
            context.hlines(0.45, gene_start, gene_end, color="#888888", linewidth=2)
            for start, end in model["exon"]:
                context.broken_barh([(start, end - start + 1)], (0.32, 0.26), facecolors="#444444")
            if len(probes):
                context.axvspan(int(probes.start.min()), int(probes.end.max()), color=color, alpha=0.2)
                for row in probes.itertuples(index=False):
                    zoom.plot([row.start, row.end], [0.7, 0.7], color=color, linewidth=6,
                              solid_capstyle="butt")
                left, right = int(probes.start.min()), int(probes.end.max())
                pad = max(100, int((right - left + 1) * 0.04))
                zoom.set_xlim(left - pad, right + pad)
                zoom.set_title(f"Selected {len(probes)} BLAST-verified probes: {right-left+1:,} bp span", loc="left")
            else:
                zoom.set_xlim(gene_start, gene_end)
                zoom.text(0.5, 0.5, "No BLAST-verified probes", ha="center", va="center",
                          transform=zoom.transAxes)
            context.set_xlim(gene_start, gene_end)
            context.set_ylim(0, 1)
            zoom.set_ylim(0, 1.4)
            context.set_yticks([])
            zoom.set_yticks([])
            context.set_title(f"{item.gene} {item.region} | {model['transcript_id']} | strand {model['strand']}", loc="left")
            context.set_ylabel("Gene")
            zoom.set_ylabel("Probes")
            for ax in (context, zoom):
                ax.xaxis.set_major_locator(MaxNLocator(nbins=5, integer=True))
                ax.xaxis.set_major_formatter(FuncFormatter(lambda value, _: f"{value:,.0f}"))
                ax.grid(axis="x", alpha=0.25)
            zoom.set_xlabel(f"{model['chrom']} genomic coordinate (mm10/GRCm38, bp)")
            fig.tight_layout()
            fig.savefig(plot_dir / f"{item.gene}_{item.region}.png", dpi=200)
            pdf.savefig(fig)
            plt.close(fig)
            bar.advance(task)
    return pdf_path
