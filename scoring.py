"""
Author: Zoe GUILBERT
Date: 2025-01-24

Description:
This script processes a PDB file to extract C3' atom data, computes distances between C3' atoms,
and calculates the estimated Gibbs free energy of the evaluated RNA conformation using linear interpolation.

Instructions:
1. Set the input file path to your PDB file.
2. Set the output folder path where you want to save the results.
3. Ensure the intervals for linear interpolation are defined correctly.
4. Run the script.

Input Table:
The input PDB file should contain atomic coordinates in the standard PDB format.
The script extracts C3' atoms and computes distances between them.

Output:
- A tabular file containing C3' atom data.
- A tabular file containing computed distances.
- The estimated Gibbs free energy printed to the console.
"""


import math
import os

# Define the input file path
input_file_path = "/home/zozo/Documents/RNA_postic/RNA_prediction/R1212TS028_1"

# Define the output folder path
output_folder_path = "/home/zozo/Documents/RNA_postic/02_results"

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
def extraction_C3_tabular_format(input_file_path, output_folder_path):
    tabular_data = []
    with open(input_file_path, 'r') as file:
        for line in file:
            parsed_line = parse_c3_prime(line)
            if parsed_line:
                tabular_data.append(parsed_line)

    # Write the results to an output file
    base_name = os.path.basename(input_file_path).replace('.pdb', '_C3_tabular.txt')
    output_file_path = os.path.join(output_folder_path, base_name)
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

# Function to compute distances between C3' atoms
def compute_distances(tabular_file_path, output_folder_path):
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
    header = ["Base Pair", "Residue i", "Residue j", "Distance (Å)"]

    # Compute distances and store them directly in tabular_data
    for i in range(len(c3_atoms)):
        for j in range(i + 4, len(c3_atoms)):  # Residues must be separated by at least 3 positions
            if c3_atoms[i][2] == c3_atoms[j][2]:  # Only consider intrachain distances
                base_pair = tuple(sorted((c3_atoms[i][0], c3_atoms[j][0])))
                if base_pair in base_pairs:
                    distance = calculate_distance(c3_atoms[i], c3_atoms[j])
                    tabular_data.append([f"{base_pair[0]}-{base_pair[1]}", c3_atoms[i][1], c3_atoms[j][1], distance])

    # Write the tabular data to a new file
    base_name = os.path.basename(tabular_file_path).replace('_C3_tabular.txt', '_distances.txt')
    output_file_path = os.path.join(output_folder_path, base_name)

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

# Main code
if __name__ == "__main__":
    # Process the PDB file and collect C3' data
    tabular_file_path = extraction_C3_tabular_format(input_file_path, output_folder_path)

    # Compute distances between C3' atoms
    distances = compute_distances(tabular_file_path, output_folder_path)

    # Define intervals for linear interpolation
    intervals = [(6, 7), (7, 8), (8, 9), (9, 10)]

    # Extract distances from the tabular data
    distance_values = [row[3] for row in distances]

    # Calculate Gibbs free energy
    gibbs_free_energy = calculate_gibbs_free_energy(distance_values, intervals)
    print(f"Estimated Gibbs Free Energy: {gibbs_free_energy}")
