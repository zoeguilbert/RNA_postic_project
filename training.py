"""
Author: Zoe GUILBERT
Date: 2025-01-24

Description:
This script processes PDB files to extract C3' atom data, computes distances between C3' atoms,
calculates the estimated energy of the evaluated RNA conformation, and plots the interaction profiles.

Instructions:
1. Set the input folder path containing PDB files.
2. Set the output folder path where you want to save the results.
3. Run the script.

Input Table:
The input PDB files should contain atomic coordinates in the standard PDB format.
The script extracts C3' atoms and computes distances between them.

Output:
- Tabular files containing C3' atom data.
- Tabular files containing computed distances.
- Tabular files containing estimated energy profiles.
- Plots of interaction profiles saved in the specified output folder.
"""

import os
import math
from collections import defaultdict
import matplotlib.pyplot as plt

# Things to Change
"""
1. Set the input folder path containing PDB files.
2. Set the output folder path where you want to save the results.
"""
RNA_pdb_folder = "/home/zozo/Documents/RNA_postic/RNA_prediction"
output_folder = "/home/zozo/Documents/RNA_postic/02_results_rna_prediction"



# Define headers for the output file
header = [
    "Residue_Name",     # Residue name (e.g., C)
    "Residue_number",   # Residue sequence number (e.g., 127)
    "Chain_ID",         # Chain ID (e.g., A)
    "X_Coordinate",     # X Coordinate (e.g., 19.212)
    "Y_Coordinate",     # Y Coordinate (e.g., 56.017)
    "Z_Coordinate"      # Z Coordinate (e.g., 58.123)
]

# Function to parse a single line for C3' atoms
def parse_c3_prime(line):
    """Parses a line to extract data if it contains a C3' atom."""
    if line.startswith("ATOM") and "C3'" in line[12:16]:
        return [
            line[17:20].strip(),               # Residue Name
            int(line[22:26].strip()),          # Residue Sequence
            line[21].strip(),                  # Chain ID
            float(line[30:38].strip()),        # X Coordinate
            float(line[38:46].strip()),        # Y Coordinate
            float(line[46:54].strip())         # Z Coordinate
        ]
    return None

# Function to process the PDB file and collect C3' data
def extraction_C3_tabular_format(input_file_path, output_folder):
    tabular_data = []
    with open(input_file_path, 'r') as file:
        for line in file:
            parsed_line = parse_c3_prime(line)
            if parsed_line:
                tabular_data.append(parsed_line)

    # Write the results to an output file
    base_name = os.path.basename(input_file_path).replace('.pdb', '_C3_tabular.txt')
    output_file_path = os.path.join(output_folder, base_name)
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')  # Write headers
        for row in tabular_data:                     # Write each parsed row
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")
    return output_file_path

# Function to calculate the distance between two points in 3D space
def calculate_distance(atom1, atom2):
    x1, y1, z1 = atom1[3], atom1[4], atom1[5]
    x2, y2, z2 = atom2[3], atom2[4], atom2[5]
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

# Function to compute distances between C3' atoms and combine results
def compute_distances_and_combine(tabular_file_path, combined_output_file_path):
    
     # Extract the base name of the file (e.g., "8olv" from "8olv_C3_tabular.txt")
    base_name = os.path.basename(tabular_file_path).replace('_C3_tabular.txt', '')
    
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
    header = ["Base_Pair", "Distance_Interval", "Observed_Frequency", "Reference_Frequency", "Log-Ratio"]

    # Compute distances and store them directly in tabular_data
    for i in range(len(c3_atoms)-4):
        for j in range(i + 4, len(c3_atoms)):  # Residues must be separated by at least 3 positions
            if c3_atoms[i][2] == c3_atoms[j][2]:  # Only consider intrachain distances
                base_pair = tuple(sorted((c3_atoms[i][0], c3_atoms[j][0])))
                if base_pair in base_pairs:
                    distance = calculate_distance(c3_atoms[i], c3_atoms[j])
                    tabular_data.append([base_name,f"{base_pair[0]}-{base_pair[1]}", c3_atoms[i][1], c3_atoms[j][1], distance])

    # Append the results to the combined output file
    with open(combined_output_file_path, 'a') as output_file:
        for row in tabular_data:
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Distances computed and combined for file: {tabular_file_path}")

# Function to calculate the estimated energy for each base pair and each distance interval
def calculate_energy(tabular_file_path, output_folder):
    # Read the tabular file and parse the data
    with open(tabular_file_path, 'r') as file:
        lines = file.readlines()

    data = [line.strip().split('\t') for line in lines[1:]] # skip header
    data = [[parts[0], parts[1], int(parts[2]),int(parts[3]), float(parts[4])] for parts in data]

    # Define distance intervals (bins)
    distance_bins = [(i, i + 2) for i in range(101)]

    # Create a dictionary to store the counts for each base pair and distance interval
    counts = defaultdict(lambda: [0] * len(distance_bins))
    total_counts = [0] * len(distance_bins)
    total_pairs = defaultdict(int)

    # Count the number of pairs for each base pair and distance interval
    for entry in data:
        rna_name, base_pair, seq_i, seq_j, distance = entry
        for i, (lower, upper) in enumerate(distance_bins):
            if lower < distance <= upper:
                counts[base_pair][i] += 1
                total_counts[i] += 1
                total_pairs[base_pair] += 1
                break
    print(counts)
    print(total_counts)
    print(total_pairs)
    # Calculate the observed frequencies
    observed_frequencies = defaultdict(lambda: [0] * len(distance_bins))
    for base_pair, count_list in counts.items():
        total_pairs_count = total_pairs[base_pair]
        for i, count in enumerate(count_list):
            if total_pairs_count > 0:
                observed_frequencies[base_pair][i] = count / total_pairs_count

    # Calculate the reference frequencies
    reference_frequencies = [0] * len(distance_bins)
    total_pairs_count = sum(total_counts)
    for i, count in enumerate(total_counts):
        if total_pairs_count > 0:
            reference_frequencies[i] = count / total_pairs_count

    # Calculate the log-ratio of the two frequencies
    log_ratios = defaultdict(lambda: [0] * len(distance_bins))
    for base_pair, freq_list in observed_frequencies.items():
        for i, freq in enumerate(freq_list):
            if freq > 0 and reference_frequencies[i] > 0:
                log_ratios[base_pair][i] = -math.log(freq / reference_frequencies[i])
            else:
                log_ratios[base_pair][i] = 10  # Set maximum scoring value to 10

    # Create tabular data for output
    tabular_data = []
    header = ["Base_Pair", "Distance_Interval", "Observed_Frequency", "Reference_Frequency", "Log-Ratio"]

    for base_pair, freq_list in observed_frequencies.items():
        for i, (lower, upper) in enumerate(distance_bins):
            tabular_data.append([
                base_pair,
                f"{lower}-{upper}",
                observed_frequencies[base_pair][i],
                reference_frequencies[i],
                log_ratios[base_pair][i]
            ])

    # Write the tabular data to a new file
    base_name = tabular_file_path.split('/')[-1].split('.')[0]
    output_file_path = f"{output_folder}/{base_name}_energy.txt"

    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for row in tabular_data:
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")
    return output_file_path


# Main function to process all PDB files in a folder
def process_all_pdb_files(folder_path, output_folder):
    # Create output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok=True)

    # Create a combined distances file
    combined_output_file_path = os.path.join(output_folder, 'combined_distances.txt')
    with open(combined_output_file_path, 'w') as output_file:
        output_file.write('\t'.join(["File_Name", "Base_Pair", "Residue_i", "Residue_j", "Distance(Å)"]) + '\n')

    # List all PDB files in the folder
    pdb_files = [f for f in os.listdir(folder_path) if f.endswith('.pdb')]

    for pdb_file in pdb_files:
        input_file_path = os.path.join(folder_path, pdb_file)

        # STEP 1: Process the PDB file and collect C3' data in tabular file
        C3_tabular_file_path = extraction_C3_tabular_format(input_file_path, output_folder)

        # STEP 2: Call the function to compute distances and combine results
        compute_distances_and_combine(C3_tabular_file_path, combined_output_file_path)

    # STEP 3: Calculate the estimated energy for each base pair and each distance interval
    energy_tab_file = calculate_energy(combined_output_file_path, output_folder)



#Main code
if __name__ == "__main__":
    # Call the main function to process all PDB files in the folder
    process_all_pdb_files(RNA_pdb_folder, output_folder)