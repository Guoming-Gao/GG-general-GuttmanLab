# output_generator.py - Focused Output Generation (FILT + HIGH_SNP only)

import pandas as pd
import os
from pathlib import Path
from rich.console import Console
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

console = Console()


class OutputGenerator:
    """Generate focused output files (FILT + HIGH_SNP only, no ALL.csv)"""

    def __init__(self, config):
        self.config = config
        self.add_rtbc = config.get("add_rtbc_barcode", True)
        self.rtbc_sequence = config.get("rtbc_sequence", "/5Phos/TGACTTGAGGAT")
        self.min_snp_coverage = config.get("min_snp_coverage_for_final", 2)
        self.generate_blast_fasta = config.get("generate_blast_fasta", True)
        self.fasta_description_format = config.get(
            "fasta_description_format", "detailed"
        )

    def generate_focused_outputs(self, filtered_probes, output_dir, file_prefix="FISH_RT_probes_FILTERED"):
        """
        Generate focused output files for high-quality probes
        """
        output_dir = Path(output_dir)
        generated_files = []

        try:
            console.print(
                f"[cyan]Generating output files for {file_prefix}...[/cyan]"
            )

            # Add RTBC barcodes to all probes
            probes_with_rtbc = self._add_rtbc_barcodes(filtered_probes)

            # Create DataFrame from filtered probes
            df_filt = pd.DataFrame(probes_with_rtbc)

            if len(df_filt) == 0:
                console.print(f"[yellow]No probes provided for {file_prefix}[/yellow]")
                return []

            # Sort by gene name, then by SNP coverage (descending), then by PNAS score
            df_filt = df_filt.sort_values(
                ["GeneName", "SNP_Count", "NbOfPNAS", "dGScore"],
                ascending=[True, False, False, False],
            )

            console.print(f"[cyan]Processed {len(df_filt)} probes[/cyan]")

            # Generate CSV probes file
            csv_file = output_dir / f"{file_prefix}.csv"
            df_filt.to_csv(csv_file, index=False)
            generated_files.append(csv_file)
            console.print(f"✅ Generated: {csv_file}")

            # Generate FASTA files
            if self.generate_blast_fasta:
                fasta_files = self._generate_focused_fasta_files(df_filt, output_dir, file_prefix)
                generated_files.extend(fasta_files)

            # Generate summary statistics (only for final selection)
            if "FINAL" in file_prefix.upper() or "FILTERED" in file_prefix.upper():
                summary_file = self._generate_focused_summary(df_filt, output_dir, file_prefix)
                if summary_file:
                    generated_files.append(summary_file)

            console.print(
                f"[green]Generated {len(generated_files)} output files[/green]"
            )
            return generated_files

        except Exception as e:
            console.print(f"[red]Error generating focused outputs: {str(e)}[/red]")
            import traceback
            console.print(traceback.format_exc())
            return []

    def _add_rtbc_barcodes(self, probe_list):
        """Add RTBC barcodes to probe sequences"""

        if not self.add_rtbc:
            for probe in probe_list:
                probe["Full_Oligo_Seq"] = probe["Probe_Seq"]
            return probe_list

        console.print(f"[cyan]Adding RTBC barcode: {self.rtbc_sequence}[/cyan]")

        for probe in probe_list:
            full_oligo_seq = self.rtbc_sequence + probe["Probe_Seq"]
            probe["Full_Oligo_Seq"] = full_oligo_seq

        return probe_list

    def _generate_focused_summary(self, df_filt, output_dir, file_prefix):
        """Generate focused summary statistics"""

        summary_file = output_dir / f"{file_prefix}_summary.txt"

        try:
            with open(summary_file, "w") as f:
                f.write(f"FISH-RT Probe Analysis Summary - {file_prefix}\n")
                f.write("=" * 50 + "\n\n")

                f.write("ANALYSIS MODE: Local Files Only (No Ensembl API)\n")
                f.write(f"Applied Filters: GC + All PNAS Rules + Dustmasker\n")
                f.write(f"Output Mode: candidates/final selection flow\n\n")

                # Overall statistics
                f.write("OVERALL STATISTICS:\n")
                f.write(f"  Probes: {len(df_filt)}\n")

                if "SNP_Count" in df_filt.columns:
                    avg_snps = df_filt["SNP_Count"].mean()
                    f.write(
                        f"  Average SNPs covered: {avg_snps:.1f}\n"
                    )

                f.write("\n")

                # Gene-level statistics
                f.write("GENE-LEVEL STATISTICS:\n")
                gene_stats = (
                    df_filt.groupby("GeneName")
                    .agg(
                        {
                            "Probe_Seq": "count",
                            "SNP_Count": "mean",
                            "dGScore": "mean",
                            "NbOfPNAS": "mean",
                        }
                    )
                    .round(2)
                )
                gene_stats.columns = [
                    "Probe_Count",
                    "Avg_SNPs",
                    "Avg_dGScore",
                    "Avg_PNAS",
                ]

                f.write(gene_stats.to_string())
                f.write("\n\n")

                # Quality metrics
                if len(df_filt) > 0:
                    f.write("QUALITY METRICS:\n")
                    f.write(f"  Average GC content: {df_filt['GC_Content'].mean():.3f}\n")
                    f.write(
                        f"  Average probe size: {df_filt['Probe_Length'].mean():.1f} nt\n"
                    )
                    f.write(f"  Average dG37 score: {df_filt['dGScore'].mean():.3f}\n")
                    f.write(
                        f"  Average PNAS score: {df_filt['NbOfPNAS'].mean():.1f}/5\n"
                    )

                    if "SNP_Count" in df_filt.columns:
                        f.write(
                            f"  Average SNPs covered: {df_filt['SNP_Count'].mean():.1f}\n"
                        )
                        f.write(
                            f"  Max SNPs covered: {df_filt['SNP_Count'].max()}\n"
                        )

                    f.write("\n")

                # SNP coverage distribution
                if "SNP_Count" in df_filt.columns:
                    f.write("SNP COVERAGE DISTRIBUTION:\n")
                    snp_counts = (
                        df_filt["SNP_Count"].value_counts().sort_index()
                    )
                    for snp_count, probe_count in snp_counts.items():
                        f.write(f"  {snp_count} SNPs: {probe_count} probes\n")
                    f.write("\n")

                # Top genes by probe count
                f.write("TOP GENES BY PROBE COUNT:\n")
                top_genes = df_filt["GeneName"].value_counts().head(10)
                for gene, count in top_genes.items():
                    f.write(f"  {gene}: {count} probes\n")

            console.print(f"✅ Generated summary: {summary_file}")
            return summary_file

        except Exception as e:
            console.print(f"[red]Error generating summary: {str(e)}[/red]")
            return None

    def _generate_focused_fasta_files(self, df_filt, output_dir, file_prefix):
        """Generate FASTA files for probe sets"""
        generated_files = []

        try:
            # Main FASTA file for the current set
            main_fasta = output_dir / f"{file_prefix}.fasta"
            self._write_fasta(df_filt, main_fasta, include_rtbc=True)
            generated_files.append(main_fasta)

            # Extra FASTA file for probe sequences only (no RTBC) - for BLAST documentation if requested
            # if self.fasta_description_format == "detailed":
            #     probes_only_fasta = output_dir / f"{file_prefix}_sequences_only.fasta"
            #     self._write_fasta(df_filt, probes_only_fasta, include_rtbc=False)
            #     generated_files.append(probes_only_fasta)

            console.print(f"✅ Generated {len(generated_files)} FASTA files")
            return generated_files

        except Exception as e:
            console.print(f"[red]Error generating FASTA files: {str(e)}[/red]")
            return []

    def _write_fasta(self, df, output_file, include_rtbc=True):
        """Write DataFrame to FASTA file"""

        records = []

        for idx, row in df.iterrows():
            # Choose sequence based on RTBC preference
            if include_rtbc and self.add_rtbc:
                sequence = row["Full_Oligo_Seq"]
                seq_type = "with_RTBC"
            else:
                sequence = row["Probe_Seq"]
                seq_type = "probe_only"

            # Create descriptive FASTA header
            gene_name = row["GeneName"]
            probe_size = row["Probe_Length"]
            snp_count = row.get("SNP_Count", 0)
            region_type = row.get("RegionType", "unknown")
            pnas_score = row.get("NbOfPNAS", 0)

            if self.fasta_description_format == "detailed":
                description = f"{gene_name} | {seq_type} | size:{probe_size}nt | SNPs:{snp_count} | region:{region_type} | PNAS:{pnas_score}/5"
            else:
                description = f"{gene_name}_{seq_type}_{idx}"

            record = SeqRecord(
                Seq(sequence), id=f"{gene_name}_probe_{idx}", description=description
            )
            records.append(record)

        SeqIO.write(records, output_file, "fasta")

