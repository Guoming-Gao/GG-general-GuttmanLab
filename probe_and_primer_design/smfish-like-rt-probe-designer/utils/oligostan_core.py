# utils/oligostan_core.py - FIXED: Correct probe sequence assignment

import numpy as np
from Bio.Seq import Seq  # ADD THIS IMPORT
from .thermodynamics import dg_calc_rna_37, dg37_score_calc
from .filters import (
    is_ok_4_pnas_filter,
    is_ok_4_gc_filter,
    is_it_ok_4_a_comp,
    is_it_ok_4_a_stack,
    is_it_ok_4_c_comp,
    is_it_ok_4_c_stack,
    is_it_ok_4_c_spec_stack,
)
from .sequence_utils import determine_region_type

def which_max_r(x):
    """Exact translation of R's WhichMax function"""
    max_val = np.max(x)
    max_indices = np.where(x == max_val)[0]
    if len(max_indices) >= 2:
        return [0, max_val]  # Return 0 for ties (R behavior)
    else:
        return [max_indices[0] + 1, max_val]  # Return 1-based index

def get_probes_from_sequence(
    seq,
    min_size_probe=26,
    max_size_probe=32,
    desired_dg=-32,
    min_score_value=0.9,
    inc_betw_prob=2,
):
    """Modified probe generation from Oligostan (same core algorithm)"""
    if isinstance(seq, list):
        seq = "".join(seq).upper()

    diff_size = max_size_probe - min_size_probe

    # Build thermodynamic matrix for all probe sizes
    the_tms_tmp = dg_calc_rna_37(seq, probe_length=max_size_probe)
    nb_of_probes = len(the_tms_tmp)

    if diff_size > 0:
        all_columns = [the_tms_tmp]  # Start with max size

        for i in range(diff_size - 1, -1, -1):
            probe_length = min_size_probe + i
            dg_values = dg_calc_rna_37(seq, probe_length=probe_length)
            min_len = min(len(dg_values), nb_of_probes)
            all_columns.insert(0, dg_values[:min_len])

        # Build matrix: rows = positions, columns = probe sizes
        max_len = max(len(col) for col in all_columns)
        the_tms_matrix = np.full((max_len, diff_size + 1), np.nan)

        for col_idx, column in enumerate(all_columns):
            the_tms_matrix[: len(column), col_idx] = column

    else:
        the_tms_matrix = np.array(the_tms_tmp).reshape(-1, 1)

    # Calculate scores
    tm_scores = np.zeros_like(the_tms_matrix)
    for col in range(the_tms_matrix.shape[1]):
        valid_mask = ~np.isnan(the_tms_matrix[:, col])
        if np.any(valid_mask):
            tm_scores[valid_mask, col] = dg37_score_calc(
                the_tms_matrix[valid_mask, col], desired_dg
            )

    # Find best probe size for each position
    best_scores = []
    for row in range(tm_scores.shape[0]):
        row_data = tm_scores[row, :]
        if not np.all(np.isnan(row_data)) and not np.all(row_data == 0):
            valid_data = row_data[~np.isnan(row_data)]
            if len(valid_data) > 0:
                result = which_max_r(valid_data)
                best_scores.append(result)
            else:
                best_scores.append([0, 0])
        else:
            best_scores.append([0, 0])

    best_scores = np.array(best_scores)

    # Convert to probe sizes
    for i in range(len(best_scores)):
        if best_scores[i, 0] != 0:  # Not a tie
            best_scores[i, 0] = best_scores[i, 0] + (min_size_probe - 1)
        else:  # Tie case - use minimum size
            best_scores[i, 0] = min_size_probe

    # Add positions
    positions = np.arange(1, len(best_scores) + 1).reshape(-1, 1)
    best_scores_with_pos = np.column_stack([best_scores, positions])

    # Filter by score threshold
    validated_scores = best_scores_with_pos[
        best_scores_with_pos[:, 1] >= min_score_value
    ]

    if len(validated_scores) == 0:
        return None

    # Apply spacing constraint
    return apply_spacing_constraint(validated_scores, seq, inc_betw_prob)

def apply_spacing_constraint(validated_scores, seq, min_spacing):
    """Apply minimum spacing between probes"""
    if len(validated_scores.shape) == 1:
        validated_scores = validated_scores.reshape(1, -1)

    validated_scores = validated_scores[np.argsort(validated_scores[:, 2])]

    the_probes = []
    pointeur = 0

    while pointeur < len(seq):
        valid_tmp = validated_scores[validated_scores[:, 2] >= pointeur]

        if len(valid_tmp) > 0:
            if len(valid_tmp.shape) == 1:
                valid_tmp = valid_tmp.reshape(1, -1)

            probe_size = int(valid_tmp[0, 0])
            score = valid_tmp[0, 1]
            position = int(valid_tmp[0, 2])

            # FIXED: Extract target sequence, then calculate probe sequence
            target_seq = seq[position - 1 : position - 1 + probe_size].upper()
            probe_seq = str(Seq(target_seq).reverse_complement())  # ACTUAL PROBE SEQUENCE

            the_probes.append(
                {
                    "size": probe_size,
                    "score": score,
                    "position": position,
                    "target_sequence": target_seq,      # NEW: target region
                    "probe_sequence": probe_seq,        # NEW: actual probe
                    "genomic_start": position,          # Will be adjusted later
                    "genomic_end": position + probe_size - 1,
                }
            )

            pointeur = position + probe_size + min_spacing
        else:
            break

    return the_probes

def design_fish_probes(gene_data, config):
    """
    Main probe design function for FISH-RT probes
    Returns list of probe dictionaries with all required information
    """
    try:
        # Extract configuration
        min_size = config["probe_length_min"]
        max_size = config["probe_length_max"]
        min_score = config["score_min"]
        spacing = config["probe_spacing_min"]
        desired_dg = config["fixed_dg37_value"]

        # Design probes using Oligostan algorithm
        raw_probes = get_probes_from_sequence(
            gene_data["sequence"],
            min_size_probe=min_size,
            max_size_probe=max_size,
            desired_dg=desired_dg,
            min_score_value=min_score,
            inc_betw_prob=spacing,
        )

        if not raw_probes:
            return []

        # Process probes and add gene context
        processed_probes = []

        for probe in raw_probes:
            # Calculate genomic coordinates (mm10)
            if gene_data["strand"] == 1:  # Forward strand
                genomic_start = gene_data["genomic_start"] + probe["position"] - 1
                genomic_end = genomic_start + probe["size"] - 1
                target_strand = "+"
                probe_strand = "-"  # Probe is complementary
            else:  # Reverse strand
                genomic_end = gene_data["genomic_end"] - probe["position"] + 1
                genomic_start = genomic_end - probe["size"] + 1
                target_strand = "-"
                probe_strand = "+"  # Probe is complementary

            # Determine if probe is in exon or intron
            # Use all_known_exons for more robust labeling if requested
            exons_to_check = gene_data.get("all_known_exons", gene_data["exon_regions"])

            region_type = determine_region_type(
                probe["position"],
                probe["position"] + probe["size"] - 1,
                exons_to_check,
                gene_data["intron_regions"],
            )

            # FIXED: Use the actual probe sequence for thermodynamic calculations
            probe_sequence = probe["probe_sequence"]  # This is the reverse complement
            target_sequence = probe["target_sequence"]  # This is the gene region

            # Calculate thermodynamic properties using PROBE sequence
            actual_dg37 = dg_calc_rna_37(
                probe_sequence, probe_length=len(probe_sequence)
            )[0]

            # Calculate quality metrics using PROBE sequence
            gc_count = probe_sequence.count("G") + probe_sequence.count("C")
            gc_percentage = gc_count / len(probe_sequence)

            # Apply quality filters using PROBE sequence
            gc_filter = is_ok_4_gc_filter(
                probe_sequence, config["gc_content_min"], config["gc_content_max"]
            )

            pnas_filter = is_ok_4_pnas_filter(
                probe_sequence, config["pnas_filter_rules"]
            )

            # Apply individual PNAS filters using PROBE sequence
            a_comp_pass = 1 if is_it_ok_4_a_comp(probe_sequence) else 0
            a_stack_pass = 1 if is_it_ok_4_a_stack(probe_sequence) else 0
            c_comp_pass = 1 if is_it_ok_4_c_comp(probe_sequence) else 0
            c_stack_pass = 1 if is_it_ok_4_c_stack(probe_sequence) else 0
            c_spec_pass = 1 if is_it_ok_4_c_spec_stack(probe_sequence) else 0

            # Calculate PNAS sum
            nb_of_pnas = (
                a_comp_pass + a_stack_pass + c_comp_pass + c_stack_pass + c_spec_pass
            )

            # Create probe record with ALL required fields
            probe_record = {
                # Gene information
                "GeneName": gene_data["gene_name"],
                "Species": "mouse",
                "Chromosome": gene_data["chromosome"],
                "RegionType": region_type,  # 'exon', 'intron', or 'exon-intron'

                # Probe sequences
                # Probe_Seq = actual oligo sequence (reverse complement of target, same as RT product)
                # Target_Seq = the mRNA region the probe hybridizes to
                "Probe_Seq": probe_sequence,
                "Target_Seq": target_sequence,
                "Probe_Length": probe["size"],

                # Probe genomic coordinates (mm10)
                "Probe_Start": genomic_start,
                "Probe_End": genomic_end,

                # Strand information
                # Target_Strand = strand of the gene/mRNA (+ or -)
                # Probe binds to mRNA, so probe sequence is complementary to mRNA
                # RT product (cDNA) has same sequence as probe
                "Target_Strand": target_strand,

                # Placeholders for SNP analysis (will be filled by snp_analyzer)
                "SNP_Count": 0,
                "SNPs_Positions": "",  # List of individual SNP positions
                "SNPs_Types": "",      # B6/Cast genotypes at each position
                "RT_Region_Start": 0,  # Start of RT extension region
                "RT_Region_End": 0,    # End of RT extension region
                "RT_Product_Seq": "",  # RT product sequence (extracted from genome)
                "Full_Oligo_Seq": "",  # Will be added by output generator (with RTBC barcode)

                # BLAST specificity (will be filled by blast_analyzer)
                "BLAST_Hits": 0,       # Number of genomic BLAST hits
                "BLAST_Identity": 0,   # % identity of best hit
                "BLAST_Unique": False, # True if exactly 1 hit (specific probe)
                "BLAST_Hit_Name": "",  # Name of unique hit target

                # Quality metrics - Oligostan inherited (kept for compatibility)
                "dG37": actual_dg37,
                "dGScore": probe["score"],
                "dGOpt": desired_dg,
                "GC_Content": gc_percentage,
                "GCFilter": 1 if gc_filter else 0,
                "PNASFilter": 1 if pnas_filter else 0,
                "aCompFilter": a_comp_pass,
                "aStackFilter": a_stack_pass,
                "cCompFilter": c_comp_pass,
                "cStackFilter": c_stack_pass,
                "cSpecStackFilter": c_spec_pass,
                "NbOfPNAS": nb_of_pnas,

                # Masking filters
                "MaskedFilter": 1,  # Default pass (dustmasker not enabled)
                "RepeatMaskerPC": 0.0,  # Default no masking
            }

            processed_probes.append(probe_record)

        return processed_probes

    except Exception as e:
        raise Exception(
            f"Error designing probes for {gene_data.get('gene_name', 'unknown')}: {str(e)}"
        )
