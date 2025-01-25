"""
Author: Zoe GUILBERT
Date: 2025-01-24

Description:
This script processes a PDB file to extract C3' atom data, computes distances between C3' atoms,
and calculates the estimated Gibbs free energy of the evaluated RNA conformation using linear interpolation.

Instructions:
1. Set the folder path containing your PDB files.
2. Set the output folder path where you want to save the results.
3. Ensure the intervals for linear interpolation are defined correctly.
4. Run the script.

Input Table:
The input PDB files should contain atomic coordinates in the standard PDB format.
The script extracts C3' atoms and computes distances between them.

Output:
- A tabular file containing C3' atom data for each pdb file.
- A tabular file containing computed distances for each pdb file.
- The estimated Gibbs free energy for each PDB file saved in a single output file.
"""

import os
from training import parse_c3_prime, extraction_C3_tabular_format, calculate_distance
import math

# Define the input file path
input_folder_path = "/home/zozo/Documents/RNA_postic/RNA_prediction/"

# Define the output folder path
output_folder_path = "/home/zozo/Documents/RNA_postic/02_results_rna_prediction"

# Ensure the output folder exists
os.makedirs(output_folder_path, exist_ok=True)

# Define headers for the output file
header = [
    "Residue Name",     # Residue name (e.g., C)
    "Residue number",   # Residue sequence number (e.g., 127)
    "Chain ID",         # Chain ID (e.g., A)
    "X Coordinate",     # X Coordinate (e.g., 19.212)
    "Y Coordinate",     # Y Coordinate (e.g., 56.017)
    "Z Coordinate"      # Z Coordinate (e.g., 58.123)
]

# Function to compute distances between C3' atoms and save results in a new file
def compute_distances(tabular_file_path, output_folder_path):
    # Extract the base name of the file (e.g., "8olv" from "8olv_C3_tabular.txt")
    base_name = os.path.basename(tabular_file_path).replace('_C3_tabular.txt', '')

    # Define the output file path
    output_file_path = os.path.join(output_folder_path, f"{base_name}_distances.txt")

    # Read the tabular file and parse the data
    with open(tabular_file_path, 'r') as file:
        lines = file.readlines()

    c3_atoms = [line.strip().split('\t') for line in lines[1:]]
    c3_atoms = [[parts[0], int(parts[1]), parts[2], float(parts[3]), float(parts[4]), float(parts[5])] for parts in c3_atoms]

    # Define valid base pairs
    base_pairs = {
        ("A", "A"), ("A", "U"), ("A", "C"), ("A", "G"),
        ("C", "C"), ("C", "G"), ("C", "U"),
        ("G", "G"), ("G", "U"),
        ("U", "U")
    }

    # Create tabular data for output
    tabular_data = []
    header = ["Base_Pair", "Residue_i", "Residue_j", "Distance (Å)"]

    # Compute distances and store them directly in tabular_data
    for i in range(len(c3_atoms) - 4):
        for j in range(i + 4, len(c3_atoms)):  # Residues must be separated by at least 3 positions
            if c3_atoms[i][2] == c3_atoms[j][2]:  # Only consider intrachain distances
                base_pair = tuple(sorted((c3_atoms[i][0], c3_atoms[j][0])))
                if base_pair in base_pairs:
                    distance = calculate_distance(c3_atoms[i], c3_atoms[j])
                    tabular_data.append([f"{base_pair[0]}-{base_pair[1]}", c3_atoms[i][1], c3_atoms[j][1], distance])

    # Write the tabular data to the new output file
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for row in tabular_data:
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")
    return tabular_data


# Function to compute the scoring value using linear interpolation
def linear_interpolation(distance, intervals):
    """
    Compute the scoring value using linear interpolation.
    """
    for interval in intervals:
        if interval[0] <= distance <= interval[1]:
            mid_point_left = (interval[0] + interval[1]) / 2 - 0.5
            mid_point_right = (interval[0] + interval[1]) / 2 + 0.5
            if distance <= mid_point_left:
                # Linear interpolation between interval[0] and mid_point_left
                return (distance - interval[0]) / (mid_point_left - interval[0])
            else:
                # Linear interpolation between mid_point_left and interval[1]
                return (interval[1] - distance) / (interval[1] - mid_point_left)
    return 0

# Function to calculate the estimated Gibbs free energy
def calculate_gibbs_free_energy(distances, intervals):
    """
    Calculate the estimated Gibbs free energy of the evaluated RNA conformation.
    """
    scores = [linear_interpolation(dist, intervals) for dist in distances]
    gibbs_free_energy = sum(scores)
    return gibbs_free_energy

# Function to process all PDB files in the input folder
def process_all_files(input_folder_path, output_folder_path, output_file):
    # Get a list of all PDB files in the input folder
    pdb_files = [f for f in os.listdir(input_folder_path) if f.endswith('.pdb')]

    with open(output_file, 'w') as outfile:
        for pdb_file in pdb_files:
            pdb_file_path = os.path.join(input_folder_path, pdb_file)

            # Process the PDB file and collect C3' data
            tabular_file_path = extraction_C3_tabular_format(pdb_file_path, output_folder_path)

            # Compute distances between C3' atoms
            distances = compute_distances(tabular_file_path, output_folder_path)

            # Define intervals for linear interpolation
            intervals = [(i, i + 1) for i in range(20)]

            # Extract distances from the tabular data
            distance_values = [row[3] for row in distances]

            # Calculate Gibbs free energy
            gibbs_free_energy = calculate_gibbs_free_energy(distance_values, intervals)
            print(f"Estimated Gibbs Free Energy for {pdb_file}: {gibbs_free_energy}")

            # Write the result to the output file
            outfile.write(f"{pdb_file}: {gibbs_free_energy}\n")

if __name__ == "__main__":
    output_file = os.path.join(output_folder_path, 'gibbs_free_energy_results.txt')
    process_all_files(input_folder_path, output_folder_path, output_file)