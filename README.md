# RNA_postic_project
Creation of an objective function for the RNA folding problem guillaume.postic@univ-evry.fr

# README

## Overview
This repository contains Python scripts designed for processing RNA PDB files and analysing interaction profiles between nucleotide pairs. The primary tasks include extracting C3' atom data, calculating distances, estimating energy profiles, and visualising interaction profiles.

## Scripts

### 1. **RNA PDB Processing and Analysis**
**File:** `training.py`

#### Description
This script processes RNA PDB files to:
- Extract C3' atom data and store it in a tabular format.
- Compute distances between C3' atoms.
- Estimate the energy profiles of RNA conformations based on distance intervals.

#### Input
- **PDB Files:** Atomic coordinates in the standard PDB format.

#### Output
- Tabular files:
  - Extracted C3' atom data.
  - Computed distances between nucleotide pairs.
  - Estimated energy profiles for nucleotide interactions.

#### Usage
1. Set the paths for the input PDB folder and output results folder:
   ```python
   RNA_pdb_folder = "/path/to/your/RNA_pdb_files"
   output_folder = "/path/to/your/output_folder"
   ```
2. Run the script:
   ```bash
   python rna_processing.py
   ```
3. Output files and plots will be saved in the specified output folder.

---

### 2. **Interaction Profiles Plotting**
**File:** `plot_interaction_profil.py`

#### Description
This script visualises the interaction profiles of nucleotide pairs based on a tabular file of interaction scores.

#### Input
- **Tabular File:**
  - Base Pair: The nucleotide pair (e.g., "A-U").
  - Distance Interval: The distance interval (e.g., "6-7").
  - Log-Ratio: The log-ratio score for the interaction.

#### Output
- Plots of interaction profiles for each base pair, saved as PNG files.

#### Usage
1. Set the paths for the input tabular file and output folder:
   ```python
   tabular_file_path = "/path/to/combined_distances_energy.txt"
   output_folder = "/path/to/interaction_profiles"
   ```
2. Run the script:
   ```bash
   python plot_interaction_profiles.py
   ```
3. Plots will be saved in the specified output folder.

---

### 3. **Scoring a predicted RNA with Gibbs free energy calculation**
**File:** `scoring.py`

#### Description
This script processes a PDB file to:
- Extract C3' atom data.
- Compute distances between C3' atoms.
- Calculate the estimated Gibbs free energy of the evaluated RNA conformation using linear interpolation.

#### Input
- **PDB File:** Atomic coordinates in the standard PDB format.

#### Output
- Tabular files:
  - C3' atom data.
  - Computed distances between nucleotide pairs.
- The estimated Gibbs free energy printed to the console.

#### Usage
1. Set the input file path and output folder:
   ```python
   input_file_path = "/path/to/your/input_pdb_file"
   output_folder_path = "/path/to/your/output_folder"
   ```
2. Run the script:
   ```bash
   python distance_energy_calculation.py
   ```
3. Output files will be saved in the specified output folder, and the estimated Gibbs free energy will be printed in the console.

---

### 4. **Nucleotide Interaction Analysis**
**File:** `Distance_comparison_big_vs_small.py`

#### Description
This script loads tabular data to:
- Categorise nucleotide pairs.
- Sort distance intervals.
- Group data by category of nucleotide pair (big-big, small-small and mixed, with big nucleotides C and G, and smaller ones A and T).
- Calculate mean log-ratios and observed frequencies.
- Visualise the results through plots.

#### Input
- **Tabular File:**
  - Base Pair: The nucleotide pair (e.g., "A-U").
  - Distance Interval: The distance interval (e.g., "6-7").
  - Observed Frequency: The observed frequency for the interaction.
  - Log-Ratio: The log-ratio score for the interaction.

#### Output
- Plots of mean log-ratios and observed frequencies for each nucleotide pair category.

#### Usage
1. Set the input file path to your tabular data file:
   ```python
   file_path = "/path/to/your/tabular_data_file.txt"
   ```
2. Run the script:
   ```bash
   python interaction_analysis.py
   ```
3. Plots will be generated and displayed or saved as needed.

## Dependencies
- Python 3.6 or higher
- Required libraries:
  - `os`
  - `math`
  - `collections`
  - `matplotlib`

Install dependencies using:
```bash
pip install matplotlib
```

## Author
**Zoe Guilbert**  
**Date:** 2025-01-24

