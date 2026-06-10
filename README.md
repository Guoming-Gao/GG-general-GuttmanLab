# GG-general-GuttmanLab Scripts

General-purpose scripts, notebooks, and small analysis packages used in the
Guttman Lab. The repository is organized by analysis domain first, with
pipeline-specific documentation kept next to the code that runs that workflow.

## Repository Layout

```text
GG-general-GuttmanLab/
|-- image_analysis/
|   |-- channel_registration/
|   |-- format_converter/
|   |-- segmentation/
|   |-- split_tif/
|   `-- pipelines/
|-- probe_and_primer_design/
|   `-- smfish-like-rt-probe-designer/
|-- sequencing_analysis/
|   |-- SPRITE/
|   |-- SWIFTseq/
|   |-- allelic_specific_amplicon_nanoporeseq/
|   |-- allelic_specific_amplicon_shortreadseq/
|   `-- focused_SWIFTseq/
|-- simulations/
`-- utilities/
```

## Main Areas

### `image_analysis/`

Image-processing tools and microscopy analysis workflows.

- Top-level `cli_*` and `gui_*` scripts provide reusable image utilities such
  as bandpass filtering, TIFF montage generation, dual-channel composites, and
  TIFF subtraction.
- `channel_registration/` contains scripts for generating and applying channel
  alignment matrices.
- `format_converter/` contains CZI, ND2, OME-TIFF, and TIFF conversion helpers,
  including MIP export scripts.
- `segmentation/` contains Cellpose/Cellpose-SAM GUI and notebook tools for
  nuclei or cell segmentation and result review.
- `split_tif/` contains channel-splitting helpers for activation- and
  camera-based acquisition formats.
- `pipelines/` contains workflow-specific notebooks and packages:
  - `BulkFluo_RDF/`: per-hub radial distribution function analysis for paired
    SPEN/H3K27ac SACD images. See
    `image_analysis/pipelines/BulkFluo_RDF/README.md`.
  - `SACDpy/`: installable `src/`-layout Python package for SACD
    reconstruction, with tests and validation reports. See
    `image_analysis/pipelines/SACDpy/README.md`.
  - `SMLM_PCF-FFT/`: FFT-accelerated SMLM pair correlation function analysis.
    See `image_analysis/pipelines/SMLM_PCF-FFT/README.md`.
  - `SMLM_PCF-old/`: older SMLM PCF notebooks kept for reference.
  - `SPT/`: single-particle tracking processing and plotting notebooks/scripts.
  - `smFISHquant/`: modular smFISH image-analysis notebooks. See
    `image_analysis/pipelines/smFISHquant/README.md`.
  - Other `nb_*` notebooks cover two-color statistics, montage generation,
    cross-section plotting, and related microscopy summaries.

### `sequencing_analysis/`

Sequencing-analysis workflows and helper tools.

- `allelic_specific_amplicon_nanoporeseq/`: Nanopore amplicon pipeline for
  whole-genome alignment followed by SNP-based B6/Cast allelic quantification.
  See `sequencing_analysis/allelic_specific_amplicon_nanoporeseq/README.md`.
- `allelic_specific_amplicon_shortreadseq/`: short-read amplicon workflow for
  SNP-based B6/Cast allelic quantification. See
  `sequencing_analysis/allelic_specific_amplicon_shortreadseq/README.md`.
- `focused_SWIFTseq/`: focused-SWIFTseq on-target rate, NRF/library complexity,
  and probe-specificity analysis. See
  `sequencing_analysis/focused_SWIFTseq/README.md`.
- `SWIFTseq/`: legacy and exploratory SWIFT-seq utilities, probe-design Java
  code, BLAST/annotation GUIs, and plotting notebooks.
- `SPRITE/`: SPRITE-related notebook work.

### `probe_and_primer_design/`

Probe and primer design workflows.

- `smfish-like-rt-probe-designer/`: focused RT-primer design using smFISH-like
  Oligostan design rules, BLAST specificity checks, and B6/Cast SNP coverage
  analysis. See
  `probe_and_primer_design/smfish-like-rt-probe-designer/README.md`.

### `utilities/`

Small reusable utilities and plotting notebooks that do not belong to one
specific pipeline.

- qPCR and fold-change plotting notebooks.
- Bar chart, dot blot, electropherogram, and AlphaFold plotting notebooks.
- `csv_to_bed.py` for simple genomic interval conversion.
- `gui_amp_cyc_calculator.py` and the accompanying cycle-number reference image.

### `simulations/`

Standalone numerical simulation notebooks and exported reports, including
normalization/significance simulations and count simulations.

## Naming Conventions

Most files follow a lightweight prefix convention:

- `cli_*.py`: command-line scripts.
- `gui_*.py`: graphical user-interface tools.
- `nb_*.ipynb`: exploratory or reusable Jupyter notebooks.
- Numbered scripts or notebooks, such as `01_*`, `02_*`, define ordered pipeline
  steps within a workflow directory.

Pipeline directories may use their own naming where that makes the workflow
clearer. For example, `SACDpy/` follows a standard Python package `src/` layout.

## Environments and Dependencies

Dependencies vary by workflow. Check the local README in each pipeline directory
before running a pipeline end to end.

Common Python dependencies across the repository include:

- `numpy`, `pandas`, `scipy`
- `matplotlib`, `seaborn`
- `scikit-image`, `tifffile`
- `pysam`, `biopython`
- `tqdm`, `rich`

Common external tools and environments include:

- `bioinfo` conda environment for sequencing workflows.
- `smlm` conda environment for SMLM, SACD, and several microscopy workflows.
- `smFISHquant` conda environment for the smFISHquant notebooks.
- Bioinformatics tools such as `samtools`, `bowtie2`, `minimap2`, `fastp`, and
  `blastn` where required.
- Microscopy packages such as `cellpose`, `spotiflow`, and `torch` where
  required.

## Working With the Repository

Use the top-level folders to find the relevant domain, then start from the
README inside the specific workflow directory when one exists. Scripts and
notebooks are generally intended to be run from their own workflow directory
unless the local documentation says otherwise.
