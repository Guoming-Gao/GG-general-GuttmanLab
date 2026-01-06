# GG-general-GuttmanLab Scripts

This repository contains general-purpose analysis scripts and notebooks developed for the Guttman Lab. The codebase has been organized into a logical structure with unified, versatile notebooks for common analysis tasks.

## Repository Structure

The repository is organized as follows:

- `analysis_pipelines/`: Complete analysis pipelines and automated processing scripts.
- `image_analysis/`: Tools for image processing, segmentation, and quantification.
    - `general_tools/`: Core image analysis tools (e.g., Two-Color Statistics).
    - `SMLM_PCF/`: Pair Correlation Function (PCF) analysis for Single Molecule Localization Microscopy.
- `visualization/`: Plotting and data visualization scripts.
- `simulations/`: Scripts for numerical simulations and modeling.
- `utilities/`: Miscellaneous helper scripts and common libraries.

## Key Notebooks

### Unified Analysis Notebooks

We have consolidated multiple task-specific notebooks into single, unified versions with configurable modes:

- **Two-Color Statistics**: `image_analysis/general_tools/nb_TwoColorStats_Montage.ipynb`
    - Supports both automated (Cellpose/intensity-based) and manual (ROI-based) segmentation.
    - Calculates distance-based statistics and colocalization metrics for two channels.
- **SMLM PCF Analysis**: `image_analysis/SMLM_PCF/nb_PCF_Unified.ipynb`
    - Combines automated cell/blob detection and manual ROI handling.
    - Features a highly optimized, vectorized calculation engine for Pair Correlation Function.
- **Bar Chart Plotting**: `visualization/nb_Plot_BarChart.ipynb`
    - Unified script for generating both vertical and horizontal bar charts with customizable aesthetics.

## Naming Conventions

All scripts and notebooks follow a standard naming convention:
`[type]_[description].[ext]`

- `cli_`: Command-line Interface script.
- `gui_`: Graphical User Interface application.
- `lib_`: Library or module containing helper functions.
- `nb_`: Jupyter Notebook.

## Requirements

The scripts in this repository primarily rely on standard scientific Python libraries:
- `numpy`, `pandas`, `scipy`, `matplotlib`, `seaborn`
- `shapely` (for geometric operations)
- `tifffile` (for image I/O)
- `tqdm` (for progress tracking)

Specific notebooks may require additional packages like `cellpose` or `scikit-image`.
