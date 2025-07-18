#!/usr/bin/env python3
"""
Oligostan: smiFISH Probe Design Tool (Python Version)
====================================================

A Python implementation of the Oligostan smiFISH probe design tool.
Designs optimal DNA/RNA probes for fluorescent in situ hybridization
experiments with thermodynamic optimization and quality control filtering.

Key Features:
- Interactive file selection dialog
- Thermodynamic optimization (ΔG37 calculations)
- Multiple quality filters (GC content, PNAS composition rules)
- Repeat masking using dustmasker (replaces RepeatMasker)
- Support for Ensembl format sequences
- Automated FLAP sequence addition for different fluorescent channels

Author: Converted from R to Python
Date: July 2025
"""

import os
import sys
import subprocess
import tempfile
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction


import argparse
from pathlib import Path
import logging
from typing import List, Dict, Tuple, Optional
import warnings

warnings.filterwarnings("ignore")

# Import tkinter for file dialog
import tkinter as tk
from tkinter import filedialog, messagebox
from tkinter import ttk

# Configure logging
logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)


def GC(sequence):
    return gc_fraction(sequence) * 100


class FileSelector:
    """GUI file selector for FASTA input"""

    def __init__(self):
        self.selected_file = None
        self.root = None

    def select_fasta_file(self):
        """Open file dialog to select FASTA file"""
        # Create a root window and hide it
        self.root = tk.Tk()
        self.root.withdraw()

        # Configure file dialog
        filetypes = [("FASTA files", "*.fasta *.fa *.fas *.seq"), ("All files", "*.*")]

        # Open file dialog
        self.selected_file = filedialog.askopenfilename(
            title="Select FASTA file for Oligostan probe design",
            filetypes=filetypes,
            initialdir=os.getcwd(),
        )

        # Close the root window
        self.root.destroy()

        if not self.selected_file:
            logger.info("No file selected. Exiting...")
            return None

        # Validate file exists
        if not os.path.exists(self.selected_file):
            messagebox.showerror(
                "Error", f"Selected file does not exist: {self.selected_file}"
            )
            return None

        # Validate file is readable
        try:
            with open(self.selected_file, "r") as f:
                first_line = f.readline().strip()
                if not first_line.startswith(">"):
                    messagebox.showwarning(
                        "Warning",
                        "Selected file may not be a valid FASTA file. "
                        "Please ensure it starts with '>' character.",
                    )
        except Exception as e:
            messagebox.showerror("Error", f"Cannot read selected file: {e}")
            return None

        logger.info(f"Selected FASTA file: {self.selected_file}")
        return self.selected_file

    def get_file_info(self):
        """Get filename and directory path from selected file"""
        if not self.selected_file:
            return None, None

        file_path = Path(self.selected_file)
        filename = file_path.name
        directory = str(file_path.parent)

        return filename, directory


class OligostanConfig:
    """Configuration class containing all Oligostan parameters"""

    def __init__(self, interactive=True):
        # File selection
        if interactive:
            self._setup_file_selection()
        else:
            # Default values for non-interactive mode
            self.fasta_filename = "mouseXist_mature.fa"
            self.fasta_path = "/Users/gmgao/Dropbox/Caltech_PostDoc_GuttmanLab/constructs_and_smiFISH/smiFISH-mXist"

        # =================== MAIN SETTINGS ===================
        # Ensembl format support
        self.ensembl_format = True

        # Thermodynamic parameters
        self.dg37_min = -36.0  # Minimum ΔG37 value to test
        self.dg37_max = -28.0  # Maximum ΔG37 value to test
        self.dg37_step = 0.5  # Step size for ΔG37 optimization
        self.dg37_desired = -32.0  # Specific ΔG37 target (if set, overrides range)
        self.score_min = 0.9  # Minimum score threshold for probes

        # Probe length parameters
        self.probe_length_min = 26  # Minimum probe length
        self.probe_length_max = 32  # Maximum probe length
        self.min_distance_between_probes = 2  # Minimum nucleotides between probes

        # Filtering parameters
        self.gc_filter = True  # Enable GC content filtering
        self.gc_min = 0.4  # Minimum GC content (40%)
        self.gc_max = 0.6  # Maximum GC content (60%)

        self.pnas_filter = False  # Enable PNAS composition rules
        self.pnas_filter_options = [1, 2, 4]  # Which PNAS rules to apply

        self.masked_filter = (
            False  # Enable repeat masking (set to True to use dustmasker)
        )
        self.max_masked_percent = 0.1  # Maximum allowed masked percentage (10%)

        # Output parameters
        self.min_probes_per_transcript = 0  # Minimum probes required per transcript

        # FLAP sequences for different fluorescent channels
        self.flap_x_seq = "CCTCCTAAGTTTCGAGCTGGACTCAGTG"
        self.flap_y_seq = "TTACACTCGGACCTCGTCGACATGCATT"
        self.flap_z_seq = "CCAGCTTCTAGCATCCATGCCCTATAAG"

        # Salt concentration for thermodynamic calculations
        self.salt_concentration = 0.115

        # Generate ΔG37 sequence for optimization
        self.dg37_sequence = list(
            np.arange(self.dg37_min, self.dg37_max + self.dg37_step, self.dg37_step)
        )

    def _setup_file_selection(self):
        """Setup file selection using GUI dialog"""
        print("=" * 60)
        print("OLIGOSTAN: smiFISH Probe Design Tool")
        print("=" * 60)
        print("Please select your FASTA file using the file dialog...")

        # Create file selector and get file
        selector = FileSelector()
        selected_file = selector.select_fasta_file()

        if not selected_file:
            print("No file selected. Exiting...")
            sys.exit(1)

        # Get filename and directory
        filename, directory = selector.get_file_info()

        if not filename or not directory:
            print("Error getting file information. Exiting...")
            sys.exit(1)

        # Set the configuration
        self.fasta_filename = filename
        self.fasta_path = directory

        print(f"Selected file: {filename}")
        print(f"Directory: {directory}")
        print("=" * 60)
        print()

    def set_file_manually(self, fasta_file_path):
        """Set FASTA file manually (for command-line use)"""
        file_path = Path(fasta_file_path)
        self.fasta_filename = file_path.name
        self.fasta_path = str(file_path.parent)

    def display_settings(self):
        """Display current configuration settings"""
        print("\nCurrent Oligostan Settings:")
        print("-" * 40)
        print(f"FASTA file: {self.fasta_filename}")
        print(f"FASTA path: {self.fasta_path}")
        print(f"Ensembl format: {self.ensembl_format}")
        print(
            f"ΔG37 range: {self.dg37_min} to {self.dg37_max} (step: {self.dg37_step})"
        )
        print(f"Probe length: {self.probe_length_min}-{self.probe_length_max} nt")
        print(f"GC filter: {self.gc_filter} ({self.gc_min}-{self.gc_max})")
        print(f"PNAS filter: {self.pnas_filter} (rules: {self.pnas_filter_options})")
        print(f"Repeat masking: {self.masked_filter}")
        print(f"Min probes per transcript: {self.min_probes_per_transcript}")
        print("-" * 40)


class DustMasker:
    """Dustmasker integration for repeat masking (replaces RepeatMasker)"""

    def __init__(self, max_masked_percent=0.1):
        self.max_masked_percent = max_masked_percent
        self.dustmasker_available = self._check_dustmasker()

    def _check_dustmasker(self):
        """Check if dustmasker is available"""
        try:
            subprocess.run(["dustmasker", "-help"], capture_output=True, check=False)
            return True
        except FileNotFoundError:
            return False

    def mask_sequences(
        self, sequences: List[str], sequence_names: List[str]
    ) -> List[float]:
        """
        Mask sequences using dustmasker and return masking percentages

        Args:
            sequences: List of DNA/RNA sequences
            sequence_names: List of sequence identifiers

        Returns:
            List of masking percentages for each sequence
        """
        if not self.dustmasker_available:
            raise RuntimeError(
                "dustmasker not found. Please install BLAST+ suite:\n"
                "conda install -c bioconda blast"
            )

        # Create temporary input file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as input_file:
            input_path = input_file.name
            for i, seq in enumerate(sequences):
                input_file.write(f">{sequence_names[i]}\n{seq}\n")

        try:
            # Run dustmasker
            output_path = input_path + ".masked"
            cmd = [
                "dustmasker",
                "-in",
                input_path,
                "-infmt",
                "fasta",
                "-out",
                output_path,
                "-outfmt",
                "fasta",
            ]

            subprocess.run(cmd, capture_output=True, text=True, check=True)

            # Calculate masking percentages
            masked_percentages = self._calculate_masked_percentages(
                output_path, len(sequences)
            )

            return masked_percentages

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"dustmasker failed: {e.stderr}")

        finally:
            # Clean up temporary files
            for file_path in [input_path, output_path]:
                if os.path.exists(file_path):
                    os.unlink(file_path)

    def _calculate_masked_percentages(
        self, masked_fasta_path: str, expected_count: int
    ) -> List[float]:
        """Calculate percentage of masked nucleotides in each sequence"""
        masked_percentages = []

        if not os.path.exists(masked_fasta_path):
            return [0.0] * expected_count

        for record in SeqIO.parse(masked_fasta_path, "fasta"):
            sequence = str(record.seq)
            total_length = len(sequence)

            # Count lowercase letters (masked regions) and 'N's
            masked_count = sum(1 for char in sequence if char.islower() or char == "N")

            # Calculate percentage
            masked_percentage = (
                (masked_count / total_length) * 100 if total_length > 0 else 0
            )
            masked_percentages.append(masked_percentage)

        return masked_percentages


class ThermodynamicCalculator:
    """Handles thermodynamic calculations for probe design"""

    def __init__(self, salt_concentration=0.115):
        self.salt_concentration = salt_concentration

        # Nearest-neighbor thermodynamic parameters (ΔG37 values in kcal/mol)
        self.nn_parameters = {
            "AA": -0.2,
            "AC": -1.5,
            "AG": -0.9,
            "AT": -1.0,
            "CA": -1.0,
            "CC": -2.2,
            "CG": -1.2,
            "CT": -1.4,
            "GA": -0.8,
            "GC": -2.4,
            "GG": -1.5,
            "GT": -1.0,
            "TA": -0.3,
            "TC": -1.4,
            "TG": -1.0,
            "TT": -0.4,
        }

    def calculate_dg37(self, sequence: str, probe_length: int) -> List[float]:
        """Calculate ΔG37 values with debugging"""
        sequence = sequence.upper()
        dg_values = []

        # Add debugging for first few probes
        for i in range(min(10, len(sequence) - probe_length + 1)):
            probe_seq = sequence[i : i + probe_length]
            dg_sum = 0

            # Sum nearest-neighbor contributions
            for j in range(len(probe_seq) - 1):
                dinucleotide = probe_seq[j : j + 2]
                if dinucleotide in self.nn_parameters:
                    dg_sum += self.nn_parameters[dinucleotide]

            # Salt concentration correction
            salt_correction = (np.log(self.salt_concentration) * -0.175) - 0.2
            dg_final = dg_sum - salt_correction

            dg_values.append(dg_final)

            # Debug first few probes
            if i < 5:
                print(f"Probe {i+1}: {probe_seq}")
                print(f"  Raw ΔG: {dg_sum:.3f}")
                print(f"  Salt correction: {salt_correction:.3f}")
                print(f"  Final ΔG37: {dg_final:.3f}")

        # Continue for rest of sequence without debugging
        for i in range(10, len(sequence) - probe_length + 1):
            probe_seq = sequence[i : i + probe_length]
            dg_sum = 0

            for j in range(len(probe_seq) - 1):
                dinucleotide = probe_seq[j : j + 2]
                if dinucleotide in self.nn_parameters:
                    dg_sum += self.nn_parameters[dinucleotide]

            salt_correction = (np.log(self.salt_concentration) * -0.175) - 0.2
            dg_final = dg_sum - salt_correction
            dg_values.append(dg_final)

        return dg_values

    def calculate_score(
        self, dg37_values: List[float], desired_dg: float = -33.0
    ) -> List[float]:
        """Calculate probe scores with debugging"""
        scores = []
        print(f"Target ΔG37: {desired_dg}")

        for i, dg in enumerate(dg37_values[:10]):  # Debug first 10
            deviation = abs(dg - desired_dg)
            score = (-0.1 * deviation) + 1
            final_score = max(0, score)
            scores.append(final_score)

            print(
                f"Probe {i+1}: ΔG37={dg:.3f}, deviation={deviation:.3f}, score={final_score:.3f}"
            )

        # Continue for rest without debugging
        for dg in dg37_values[10:]:
            deviation = abs(dg - desired_dg)
            score = (-0.1 * deviation) + 1
            scores.append(max(0, score))

        return scores

    # def debug_thermodynamic_calculation(self, sequence: str, probe_length: int = 26):
    #     """Debug function to compare with R version"""
    #     sequence = sequence.upper()

    #     print(f"=== Debugging Thermodynamic Calculation ===")
    #     print(f"Sequence: {sequence}")
    #     print(f"Probe length: {probe_length}")

    #     # Calculate first few probes manually
    #     for i in range(min(5, len(sequence) - probe_length + 1)):
    #         probe_seq = sequence[i:i + probe_length]
    #         print(f"\nProbe {i+1}: {probe_seq}")

    #         # Calculate dinucleotides
    #         dg_sum = 0
    #         dinucleotides = []
    #         for j in range(len(probe_seq) - 1):
    #             dinuc = probe_seq[j:j + 2]
    #             dg_val = self.nn_parameters.get(dinuc, 0)
    #             dg_sum += dg_val
    #             dinucleotides.append(f"{dinuc}:{dg_val}")

    #         print(f"  Dinucleotides: {', '.join(dinucleotides)}")
    #         print(f"  Raw dG sum: {dg_sum:.3f}")

    #         # Salt correction
    #         salt_correction = (np.log(self.salt_concentration) * -0.175) - 0.2
    #         dg_final = dg_sum - salt_correction

    #         print(f"  Salt correction: {salt_correction:.3f}")
    #         print(f"  Final dG37: {dg_final:.3f}")

    #         # Score calculation
    #         deviation = abs(dg_final - (-32.0))
    #         score = max(0, (-0.1 * deviation) + 1)
    #         print(f"  Score (target -32.0): {score:.3f}")


class ProbeFilter:
    """Handles all probe filtering operations"""

    def __init__(self, config: OligostanConfig):
        self.config = config
        self.dustmasker = (
            DustMasker(config.max_masked_percent) if config.masked_filter else None
        )

    def filter_gc_content(self, sequences: List[str]) -> List[bool]:
        """Filter probes based on GC content"""
        if not self.config.gc_filter:
            return [True] * len(sequences)

        passed = []
        for seq in sequences:
            gc_content = GC(seq) / 100.0  # Convert percentage to fraction
            passes = self.config.gc_min <= gc_content <= self.config.gc_max
            passed.append(passes)

        return passed

    def filter_pnas_rules(self, sequences: List[str]) -> List[bool]:
        """Apply PNAS composition rules"""
        if not self.config.pnas_filter:
            return [True] * len(sequences)

        passed = []
        for seq in sequences:
            passes = self._apply_pnas_rules(seq, self.config.pnas_filter_options)
            passed.append(passes)

        return passed

    def _apply_pnas_rules(self, sequence: str, rules_to_apply: List[int]) -> bool:
        """Apply individual PNAS rules to a sequence"""
        sequence = sequence.upper()

        # Rule 1: A content < 28%
        if 1 in rules_to_apply:
            a_content = sequence.count("A") / len(sequence)
            if a_content >= 0.28:
                return False

        # Rule 2: No AAAA repeats
        if 2 in rules_to_apply:
            if "AAAA" in sequence:
                return False

        # Rule 3: C content between 22% and 28%
        if 3 in rules_to_apply:
            c_content = sequence.count("C") / len(sequence)
            if not (0.22 < c_content < 0.28):
                return False

        # Rule 4: No CCCC repeats
        if 4 in rules_to_apply:
            if "CCCC" in sequence:
                return False

        # Rule 5: No 6-nucleotide window with >50% C content
        if 5 in rules_to_apply:
            for i in range(len(sequence) - 5):
                window = sequence[i : i + 6]
                c_content_window = window.count("C") / 6
                if c_content_window > 0.5:
                    return False

        return True

    def filter_repeat_masking(
        self, sequences: List[str]
    ) -> Tuple[List[float], List[bool]]:
        """Apply repeat masking filter using dustmasker"""
        if not self.config.masked_filter:
            return [0.0] * len(sequences), [True] * len(sequences)

        sequence_names = [f"probe_{i}" for i in range(len(sequences))]

        try:
            masked_percentages = self.dustmasker.mask_sequences(
                sequences, sequence_names
            )
            max_percent = self.config.max_masked_percent * 100
            passed = [percentage < max_percent for percentage in masked_percentages]

            return masked_percentages, passed

        except RuntimeError as e:
            logger.error(f"Repeat masking failed: {e}")
            raise


class ProbeGenerator:
    """Main probe generation engine"""

    def __init__(self, config: OligostanConfig):
        self.config = config
        self.thermo_calc = ThermodynamicCalculator(config.salt_concentration)
        self.probe_filter = ProbeFilter(config)

    def generate_probes_from_sequence(
        self, sequence: str, sequence_name: str, desired_dg: float = -33.0
    ) -> pd.DataFrame:
        """
        Generate optimal probes from a single sequence (matches R version logic)
        """
        sequence = sequence.upper()

        # Calculate dG37 for the longest probe length first
        max_length = self.config.probe_length_max
        min_length = self.config.probe_length_min

        # Get dG37 values for maximum probe length
        dg37_max = self.thermo_calc.calculate_dg37(sequence, max_length)

        if not dg37_max:
            return pd.DataFrame()

        num_positions = len(dg37_max)

        # Create matrix for all probe lengths (like R version)
        dg37_matrix = np.zeros((num_positions, max_length - min_length + 1))

        # Fill matrix with dG37 values for each probe length
        for i, probe_length in enumerate(range(min_length, max_length + 1)):
            dg37_values = self.thermo_calc.calculate_dg37(sequence, probe_length)
            # Pad with zeros if needed (shorter probes have fewer positions)
            if len(dg37_values) < num_positions:
                dg37_values.extend([0] * (num_positions - len(dg37_values)))
            dg37_matrix[: len(dg37_values), i] = dg37_values[:num_positions]

        # Calculate scores for all positions and lengths
        score_matrix = np.zeros_like(dg37_matrix)
        for i in range(dg37_matrix.shape[1]):
            for j in range(dg37_matrix.shape[0]):
                if dg37_matrix[j, i] != 0:  # Only calculate for valid positions
                    deviation = abs(dg37_matrix[j, i] - desired_dg)
                    score = max(0, (-0.1 * deviation) + 1)
                    score_matrix[j, i] = score

        # Select best probe length for each position (WhichMax equivalent)
        all_probes = []
        for pos in range(num_positions):
            # Find best score for this position
            best_score_idx = np.argmax(score_matrix[pos, :])
            best_score = score_matrix[pos, best_score_idx]

            # Only keep if score is above threshold
            if best_score >= self.config.score_min:
                best_length = min_length + best_score_idx
                best_dg37 = dg37_matrix[pos, best_score_idx]

                # Extract the actual probe sequence
                probe_seq = sequence[pos : pos + best_length]

                probe_info = {
                    "sequence_name": sequence_name,
                    "probe_sequence": probe_seq,
                    "start_pos": pos + 1,  # 1-based indexing
                    "end_pos": pos + best_length,
                    "probe_length": best_length,
                    "dg37": best_dg37,
                    "score": best_score,
                    "desired_dg": desired_dg,
                }

                all_probes.append(probe_info)

        if not all_probes:
            return pd.DataFrame()

        # Convert to DataFrame and sort by position (like R version)
        probes_df = pd.DataFrame(all_probes)
        probes_df = probes_df.sort_values("start_pos")

        # Apply spacing constraint
        selected_probes = self._apply_spacing_constraint(probes_df)

        return selected_probes

    def _apply_spacing_constraint(self, probes_df: pd.DataFrame) -> pd.DataFrame:
        """Apply minimum distance constraint between probes"""
        if probes_df.empty:
            return probes_df

        selected_indices = []
        last_end_pos = 0

        for idx, row in probes_df.iterrows():
            if (
                row["start_pos"]
                >= last_end_pos + self.config.min_distance_between_probes
            ):
                selected_indices.append(idx)
                last_end_pos = row["end_pos"]

        return probes_df.loc[selected_indices].reset_index(drop=True)

    def apply_all_filters(self, probes_df: pd.DataFrame) -> pd.DataFrame:
        """Apply all filtering steps to probe DataFrame"""
        if probes_df.empty:
            return probes_df

        sequences = probes_df["probe_sequence"].tolist()

        # Apply GC content filter
        gc_passed = self.probe_filter.filter_gc_content(sequences)
        probes_df["gc_content"] = [GC(seq) for seq in sequences]
        probes_df["gc_filter"] = gc_passed

        # Apply PNAS rules filter
        pnas_passed = self.probe_filter.filter_pnas_rules(sequences)
        probes_df["pnas_filter"] = pnas_passed

        # Apply repeat masking filter
        if self.config.masked_filter:
            masked_percentages, masked_passed = self.probe_filter.filter_repeat_masking(
                sequences
            )
            probes_df["masked_percentage"] = masked_percentages
            probes_df["masked_filter"] = masked_passed
        else:
            probes_df["masked_percentage"] = 0.0
            probes_df["masked_filter"] = True

        # Calculate overall filter
        probes_df["overall_filter"] = (
            probes_df["gc_filter"]
            & probes_df["pnas_filter"]
            & probes_df["masked_filter"]
        )

        return probes_df


class OligostanPipeline:
    """Main Oligostan pipeline orchestrator"""

    def __init__(self, config: OligostanConfig):
        self.config = config
        self.probe_generator = ProbeGenerator(config)

    def load_sequences(self, fasta_path: str) -> List[Tuple[str, str]]:
        """Load sequences from FASTA file"""
        sequences = []

        try:
            for record in SeqIO.parse(fasta_path, "fasta"):
                # Reverse complement for probe perspective (matching R script)
                seq = str(record.seq.reverse_complement())
                name = record.id
                sequences.append((name, seq))

            logger.info(f"Loaded {len(sequences)} sequences from {fasta_path}")
            return sequences

        except Exception as e:
            logger.error(f"Error loading sequences: {e}")
            raise

    def optimize_dg37(self, sequences: List[Tuple[str, str]]) -> float:
        """
        Optimize ΔG37 value by testing different values and selecting the best

        Args:
            sequences: List of (name, sequence) tuples

        Returns:
            Optimal ΔG37 value
        """
        if self.config.dg37_desired is not None:
            return self.config.dg37_desired

        logger.info("Optimizing ΔG37 value...")

        dg37_results = {}

        for dg37 in self.config.dg37_sequence:
            total_probes = 0
            total_transcripts_with_min_probes = 0

            for seq_name, seq in sequences:
                probes_df = self.probe_generator.generate_probes_from_sequence(
                    seq, seq_name, dg37
                )
                probes_df = self.probe_generator.apply_all_filters(probes_df)

                passed_probes = probes_df[probes_df["overall_filter"]]
                total_probes += len(passed_probes)

                if len(passed_probes) >= self.config.min_probes_per_transcript:
                    total_transcripts_with_min_probes += 1

            dg37_results[dg37] = {
                "total_probes": total_probes,
                "transcripts_with_min_probes": total_transcripts_with_min_probes,
            }

        # Select ΔG37 with maximum transcripts having minimum probes,
        # then maximum total probes as tiebreaker
        best_dg37 = max(
            dg37_results.keys(),
            key=lambda x: (
                dg37_results[x]["transcripts_with_min_probes"],
                dg37_results[x]["total_probes"],
            ),
        )

        logger.info(f"Selected optimal ΔG37: {best_dg37}")
        return best_dg37

    def add_flap_sequences(self, probes_df: pd.DataFrame) -> pd.DataFrame:
        """Add FLAP sequences for different fluorescent channels"""
        if probes_df.empty:
            return probes_df

        probes_df = probes_df.copy()

        probes_df["flap_x"] = probes_df["probe_sequence"] + self.config.flap_x_seq
        probes_df["flap_y"] = probes_df["probe_sequence"] + self.config.flap_y_seq
        probes_df["flap_z"] = probes_df["probe_sequence"] + self.config.flap_z_seq

        return probes_df

    def run(self, output_dir: str = None) -> Dict[str, pd.DataFrame]:
        """
        Run the complete Oligostan pipeline

        Args:
            output_dir: Output directory for results

        Returns:
            Dictionary containing raw and filtered probe DataFrames
        """
        logger.info("Starting Oligostan pipeline...")

        # Display current settings
        self.config.display_settings()

        # Setup output directory - CREATE IN SAME DIRECTORY AS FASTA FILE
        if output_dir is None:
            base_name = Path(self.config.fasta_filename).stem
            # Create output directory in the same location as the FASTA file
            output_dir = os.path.join(self.config.fasta_path, f"Probes_{base_name}")

        os.makedirs(output_dir, exist_ok=True)

        # Print where results will be saved
        logger.info(f"Output directory: {output_dir}")

        # Load sequences
        fasta_path = os.path.join(self.config.fasta_path, self.config.fasta_filename)
        sequences = self.load_sequences(fasta_path)

        # Optimize ΔG37
        optimal_dg37 = self.optimize_dg37(sequences)

        # Generate probes for all sequences
        all_probes = []

        for seq_name, seq in sequences:
            logger.info(f"Processing sequence: {seq_name}")

            probes_df = self.probe_generator.generate_probes_from_sequence(
                seq, seq_name, optimal_dg37
            )
            probes_df = self.probe_generator.apply_all_filters(probes_df)

            if not probes_df.empty:
                all_probes.append(probes_df)

        # Combine all probes
        if all_probes:
            combined_probes = pd.concat(all_probes, ignore_index=True)
        else:
            combined_probes = pd.DataFrame()

        # Add FLAP sequences
        combined_probes = self.add_flap_sequences(combined_probes)

        # Filter for final results
        filtered_probes = combined_probes[combined_probes["overall_filter"]].copy()

        # Save results
        self._save_results(combined_probes, filtered_probes, output_dir)

        results = {
            "raw_probes": combined_probes,
            "filtered_probes": filtered_probes,
            "optimal_dg37": optimal_dg37,
        }

        logger.info("Oligostan pipeline completed successfully!")
        return results

    def _save_results(
        self, raw_probes: pd.DataFrame, filtered_probes: pd.DataFrame, output_dir: str
    ):
        """Save results to output files"""
        base_name = Path(self.config.fasta_filename).stem

        # Save raw probes
        raw_filename = os.path.join(output_dir, f"Probes_{base_name}_ALL.txt")
        raw_probes.to_csv(raw_filename, sep="\t", index=False)

        # Save filtered probes
        filtered_filename = os.path.join(output_dir, f"Probes_{base_name}_FILT.txt")
        if not filtered_probes.empty:
            filtered_probes.to_csv(filtered_filename, sep="\t", index=False)
        else:
            with open(filtered_filename, "w") as f:
                f.write(
                    "No probes found after filtering. Consider adjusting parameters.\n"
                )

        # Save FASTA files
        self._save_fasta(
            raw_probes,
            os.path.join(output_dir, f"Probes_{base_name}_ALL_summary.fasta"),
        )
        self._save_fasta(
            filtered_probes,
            os.path.join(output_dir, f"Probes_{base_name}_FILT_summary.fasta"),
        )

        # Save settings
        settings_file = os.path.join(output_dir, "Oligostan_Settings.txt")
        with open(settings_file, "w") as f:
            f.write("Oligostan smiFISH probe design (Python version)\n")
            f.write(f"Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"Input file: {self.config.fasta_filename}\n")
            f.write(f"Input path: {self.config.fasta_path}\n")
            f.write(f"Total raw probes: {len(raw_probes)}\n")
            f.write(f"Total filtered probes: {len(filtered_probes)}\n")

        logger.info(f"Results saved to {output_dir}")

    def _save_fasta(self, probes_df: pd.DataFrame, filename: str):
        """Save probes as FASTA file"""
        if probes_df.empty:
            return

        records = []
        for _, row in probes_df.iterrows():
            record_id = f"{row['sequence_name']} {row['start_pos']} {row['end_pos']}"
            record = SeqRecord(Seq(row["probe_sequence"]), id=record_id, description="")
            records.append(record)

        SeqIO.write(records, filename, "fasta")


def main():
    """Main function with command-line interface"""
    parser = argparse.ArgumentParser(description="Oligostan: smiFISH Probe Design Tool")
    parser.add_argument(
        "--fasta", help="Input FASTA file path (overrides GUI selection)"
    )
    parser.add_argument("--output", help="Output directory")
    parser.add_argument("--dg37", type=float, help="Specific ΔG37 value to use")
    parser.add_argument(
        "--no-masking", action="store_true", help="Disable repeat masking"
    )
    parser.add_argument(
        "--min-probes", type=int, default=0, help="Minimum probes per transcript"
    )
    parser.add_argument(
        "--no-gui",
        action="store_true",
        help="Disable GUI file selection (requires --fasta)",
    )

    args = parser.parse_args()

    # Create configuration
    if args.fasta and args.no_gui:
        # Command-line mode
        config = OligostanConfig(interactive=False)
        config.set_file_manually(args.fasta)
    elif args.fasta and not args.no_gui:
        # Use command-line file but still show settings
        config = OligostanConfig(interactive=False)
        config.set_file_manually(args.fasta)
    else:
        # Interactive mode with GUI
        config = OligostanConfig(interactive=True)

    # Apply command-line overrides
    if args.dg37:
        config.dg37_desired = args.dg37

    if args.no_masking:
        config.masked_filter = False

    config.min_probes_per_transcript = args.min_probes

    # Run pipeline
    pipeline = OligostanPipeline(config)
    results = pipeline.run(args.output)

    # Print summary
    print(f"\n{'='*60}")
    print("OLIGOSTAN RESULTS SUMMARY")
    print(f"{'='*60}")
    print(f"Input file: {config.fasta_filename}")
    print(f"Total raw probes generated: {len(results['raw_probes'])}")
    print(f"Total filtered probes: {len(results['filtered_probes'])}")
    print(f"Optimal ΔG37 value: {results['optimal_dg37']:.1f}")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
