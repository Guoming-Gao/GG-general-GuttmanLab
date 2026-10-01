# filters.py
import subprocess
import tempfile
import os
from Bio import SeqIO
from Bio.Seq import Seq

try:
    from Bio.SeqUtils import gc_fraction
except ImportError:
    from Bio.SeqUtils import GC

    def gc_fraction(seq):
        return GC(seq) / 100


def is_ok_4_pnas_filter(seq, filter_to_be_used=[1, 2, 3, 4, 5]):
    """Exact translation of isOk4PNASFilter from R"""
    results = []

    if 1 in filter_to_be_used:
        results.append(is_it_ok_4_a_comp(seq))
    if 2 in filter_to_be_used:
        results.append(is_it_ok_4_a_stack(seq))
    if 3 in filter_to_be_used:
        results.append(is_it_ok_4_c_comp(seq))
    if 4 in filter_to_be_used:
        results.append(is_it_ok_4_c_stack(seq))
    if 5 in filter_to_be_used:
        results.append(is_it_ok_4_c_spec_stack(seq))

    return all(results)


def is_it_ok_4_a_comp(seq):
    """PNAS Rule 1: Adenine content < 28%"""
    a_count = seq.upper().count("A")
    return (a_count / len(seq)) < 0.28


def is_it_ok_4_a_stack(seq):
    """PNAS Rule 2: No AAAA runs"""
    return "AAAA" not in seq.upper()


def is_it_ok_4_c_comp(seq):
    """PNAS Rule 3: Cytosine content between 22-28%"""
    c_count = seq.upper().count("C")
    c_comp = c_count / len(seq)
    return 0.22 < c_comp < 0.28


def is_it_ok_4_c_stack(seq):
    """PNAS Rule 4: No CCCC runs"""
    return "CCCC" not in seq.upper()


def is_it_ok_4_c_spec_stack(seq):
    """Match Oligostan.r Rule 5, which examines its first six 6-nt windows.

    The original R script builds a fixed 6x6 index matrix, so it does not
    scan windows beyond positions 1-6. Keep this behavior for output parity.
    """
    seq = seq.upper()
    for i in range(min(6, max(0, len(seq) - 5))):
        window = seq[i : i + 6]
        c_percent = window.count("C") / 6
        if c_percent > 0.5:
            return False
    return True


def longest_homopolymer_run(seq):
    """Length of the longest identical-base run in a probe sequence."""
    if not seq:
        return 0
    seq = seq.upper()
    longest = current = 1
    for previous, base in zip(seq, seq[1:]):
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
    return longest


def is_ok_4_homopolymer(seq, max_len=4):
    """Match the RT designer: reject A/T/C/G runs longer than max_len."""
    if max_len < 1:
        raise ValueError("max_len must be positive")
    seq = seq.upper()
    return not any(base * (max_len + 1) in seq for base in "ATCG")


def is_ok_4_gc_filter(seq, min_gc=0.4, max_gc=0.6):
    """GC content filter matching R script logic"""
    gc_content = gc_fraction(seq.upper())
    return min_gc <= gc_content <= max_gc


def dustmasker_filter(sequences, max_masked_percent=0.1):
    """
    RESTORED: dustmasker filter to replace RepeatMasker
    Returns tuple: (filter_results, masked_percentages)
    filter_results: list of booleans (True = pass, False = fail)
    masked_percentages: list of masked percentages for each sequence
    """
    if not sequences:
        return [], []

    # Create temporary input file
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as temp_input:
        for i, seq in enumerate(sequences):
            temp_input.write(f">seq_{i}\n{seq}\n")
        temp_input_path = temp_input.name

    temp_output_path = temp_input_path + ".masked"

    try:
        # Run dustmasker
        result = subprocess.run(
            [
                "dustmasker",
                "-in",
                temp_input_path,
                "-out",
                temp_output_path,
                "-outfmt",
                "fasta",
            ],
            capture_output=True,
            text=True,
        )

        # Parse masked sequences and calculate masked percentages
        masked_percentages = []
        filter_results = []

        if result.returncode == 0 and os.path.exists(temp_output_path):
            # Parse dustmasker output
            masked_sequences = []
            for record in SeqIO.parse(temp_output_path, "fasta"):
                masked_sequences.append(str(record.seq))

            # Calculate masked percentages
            for i, (original_seq, masked_seq) in enumerate(
                zip(sequences, masked_sequences)
            ):
                if len(masked_seq) > 0:
                    # Count lowercase nucleotides (masked regions)
                    total_length = len(masked_seq)
                    masked_count = sum(1 for char in masked_seq if char.islower())
                    masked_percent = (
                        masked_count / total_length if total_length > 0 else 0
                    )
                else:
                    masked_percent = 0

                masked_percentages.append(masked_percent)
                # Pass filter if masked percentage is below threshold
                filter_results.append(masked_percent <= max_masked_percent)

        else:
            # dustmasker failed - pass all sequences (like R script behavior)
            print(f"Warning: dustmasker failed (return code: {result.returncode})")
            print(f"stderr: {result.stderr}")
            masked_percentages = [0.0] * len(sequences)
            filter_results = [True] * len(sequences)

    except FileNotFoundError:
        # dustmasker not found - pass all sequences (graceful degradation)
        print("Warning: dustmasker not found. Skipping repeat masking filter.")
        masked_percentages = [0.0] * len(sequences)
        filter_results = [True] * len(sequences)

    except Exception as e:
        # Any other error - pass all sequences
        print(f"Warning: dustmasker error: {e}")
        masked_percentages = [0.0] * len(sequences)
        filter_results = [True] * len(sequences)

    finally:
        # Clean up temp files
        if os.path.exists(temp_input_path):
            os.unlink(temp_input_path)
        if os.path.exists(temp_output_path):
            os.unlink(temp_output_path)

    return filter_results, masked_percentages
