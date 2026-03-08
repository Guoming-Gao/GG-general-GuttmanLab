# smFISHquant Pipeline

This repository contains a modular pipeline for multi-color smFISH image analysis.

## Project Structure

- `01_ZStack_Compiler.ipynb`: Compiles raw ONI TIFs into Z-stacks and MIPs.
- `02_Nuclei_Segmentation.ipynb`: Segments nuclei from DAPI channel.
- `03_Spot_Detection.ipynb`: Detects spots and assigns them to nuclei.
- `04_Colocalization_Analysis.ipynb`: Performs colocalization and summarizes results.
- `smfishquant/`: Core Python package containing shared logic.
- `environment.yml`: Conda environment specification.

## Setup Instructions

1.  **Create the environment**:
    ```bash
    mamba env create -f environment.yml
    ```
    *Note: This will also install `spotiflow` from the local path `/Users/gmgao/GGscripts/spotiflow`.*

2.  **Activate the environment**:
    ```bash
    conda activate smFISHquant
    ```

3.  **Register Jupyter Kernel**:
    ```bash
    python -m ipykernel install --user --name smfishquant --display-name "Python (smFISHquant)"
    ```

4.  **Run the notebooks**:
    Open the notebooks in Jupyter Lab/Notebook and ensure you select the `smFISHquant` kernel. Each notebook contains all necessary functions and logic for its specific step, allowing for easy parameter tuning and visibility.

