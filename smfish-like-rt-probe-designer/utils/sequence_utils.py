# utils/sequence_utils.py - FIXED: Removed automatic reverse complement

from Bio import SeqIO
from Bio.Seq import Seq
import os


def read_fasta_sequences(file_path):
    """
    FIXED: Read FASTA sequences without automatic reverse complement

    Now respects the sequence orientation as saved by gene_fetcher.py
    (which should already be the correct gene sequence)
    """
    sequences = []
    # Extract base filename without extension
    base_filename = os.path.splitext(os.path.basename(file_path))[0]

    for record in SeqIO.parse(file_path, "fasta"):
        # Use base filename instead of sequence header
        seq_id = base_filename

        # FIXED: Use sequence as-is (should already be correct gene sequence)
        # No more automatic reverse complement
        sequence = str(record.seq)

        sequences.append({"id": seq_id, "name": base_filename, "sequence": sequence})

    return sequences


def determine_region_type(probe_start, probe_end, exon_regions, intron_regions):
    """
    Determine if probe is in exon, intron, or spans both.
    Handles 1-indexed inclusive coordinates.
    Returns: 'exon', 'intron', or 'exon-intron'
    """
    # Check overlap with exons
    exon_overlap = False
    for exon in exon_regions:
        # Standard range overlap check: [s1, e1] overlaps [s2, e2] if s1 <= e2 and e1 >= s2
        if probe_start <= exon["end"] and probe_end >= exon["start"]:
            exon_overlap = True
            break

    # Check overlap with introns
    intron_overlap = False
    for intron in intron_regions:
        if probe_start <= intron["end"] and probe_end >= intron["start"]:
            intron_overlap = True
            break

    # Determine region type
    if exon_overlap and intron_overlap:
        return "exon-intron"
    elif exon_overlap:
        return "exon"
    elif intron_overlap:
        return "intron"
    else:
        return "intergenic"  # Shouldn't happen for gene-based design


def create_output_directory(input_file_path):
    """Create output directory in same location as input file"""
    input_dir = os.path.dirname(input_file_path)
    base_name = os.path.splitext(os.path.basename(input_file_path))[0]
    output_dir = os.path.join(input_dir, f"Probes_{base_name}")

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return output_dir
