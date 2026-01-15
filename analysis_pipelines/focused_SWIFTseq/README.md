# Focused-SWIFTseq Analysis Pipeline

This pipeline provides a robust framework for quantifying **on-target rates** and **library complexity** (NRF) from focused-SWIFTseq datasets, specifically optimized for RT-primed fragments.

## Pipeline Overview

The analysis is structured into 6 modular Python scripts:

1.  **`step1_process_probes.py`**: Loads the FISH-RT probe list, assigns genomic strands from reference GTF (mm10), and defines the refined **Initiation Bins**.
2.  **`step2_on_target_rates.py`**: Quantifies reads mapping to the initiation bin using **Strict 5'-End Detection**.
3.  **`step3_complexity.py`**: Calculates the **Non-Redundant Fraction (NRF)** as the ratio of unique sequencing events (Deduplicated) to total events (Raw) for on-target reads.
4.  **`step4_summary.py`**: Aggregates metrics across all conditions and generates the **Master Summary Report**.
5.  **`step5_spot_checks.py`**: Generates comparative paired histograms for every probe to allow for visual audit.
6.  **`step6_filter_probes.py`**: Performs an automated **Specificity Audit** using peak detection relative to a local background.

## Key Technical Considerations

### 1. RT Initiation Directionality
A critical discovery in this refinement was the correction of the RT elongation direction:
- **`+` Strand RNA**: cDNA extends towards lower genomic coordinates. The initiation bin is placed **upstream** (`Start - 15` to `Start - 1`).
- **`-` Strand RNA**: cDNA extends towards higher genomic coordinates. The initiation bin is placed **downstream** (`End + 1` to `End + 15`).

### 2. Strict Initiation Detection
To eliminate background noise from long transcripts passing through the target region, the pipeline implements a **Strict 5' End Check**. A read is classified as "on-target" ONLY if its biological 5' end (initiation site) falls within the 1-15 nt bin.

### 3. Local Background Filtering
Probes are filtered for specificity by comparing the initiation peak height to the **Mode** of the read start counts in a **±500 bp background**. Probes are only "Accepted" if they show a clear localized signal above this baseline.

### 4. Consolidated Result Management
To ensure data portability, all results (intermediate CSVs, reports, and plots) are stored in a centralized `results/` folder within the dataset directory:
`.../20251214_vs_1218-focusedSwift-RTprimerconc_UMIsplintconc/results/`

## Requirements
- `pysam`
- `pandas`
- `matplotlib`
- `seaborn`
- `scipy`
- `bioinfo` conda environment
