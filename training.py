# bla



## Part I : Extracting C3' coordinates and informations


### STEP 1 : Process the PDB file and collect C3' data

# Define the input file path
input_file_path = "/home/zozo/Documents/RNA_postic/RNA_pdb/1gid.pdb"

# Define headers for the output file
header = [
    "Residue Name",     # Residue name (e.g., C)
    "Residue number", # Residue sequence number (e.g., 127)
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

# Function Process the PDB file and collect C3' data
def extraction_C3_tabular_format(input_file_path):
    tabular_data = []
    with open(input_file_path, 'r') as file:
        for line in file:
            parsed_line = parse_c3_prime(line)
            if parsed_line:
                print(parsed_line)
                tabular_data.append(parsed_line)

    # Write the results to an output file
    output_file_path = input_file_path.split('/')[-1].replace('.pdb', '_C3_tabular.txt')
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')  # Write headers
        for row in tabular_data:                     # Write each parsed row
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")
    return output_file_path

### STEP 2 : Calculation of distances

import math

# Function to calculate the distance between two points in 3D space
def calculate_distance(atom1, atom2):
    x1, y1, z1 = atom1[3], atom1[4], atom1[5]
    x2, y2, z2 = atom2[3], atom2[4], atom2[5]
    return math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2 + (z2 - z1) ** 2)

# Function to compute distances between C3' atoms
def compute_distances(tabular_file_path):
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
    base_name = tabular_file_path.split('/')[-1].split('.')[0]
    output_file_path = f"{base_name}_distances.txt"

    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for row in tabular_data:
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")
    return output_file_path
 

# STEP 3 : Calculate the estimated enregy for each base pairs and each distances intervals

import math
from collections import defaultdict

# Function to calculate the estimated energy for each base pair and each distance interval
def calculate_energy(tabular_file_path):
    # Read the tabular file and parse the data
    with open(tabular_file_path, 'r') as file:
        lines = file.readlines()

    data = [line.strip().split('\t') for line in lines[1:]]
    data = [[parts[0], int(parts[1]), int(parts[2]), float(parts[3])] for parts in data]

    # Define distance intervals (bins)
    distance_bins = [(i, i + 1) for i in range(21)]

    # Create a dictionary to store the counts for each base pair and distance interval
    counts = defaultdict(lambda: [0] * len(distance_bins))
    total_counts = [0] * len(distance_bins)
    total_pairs = defaultdict(int)

    # Count the number of pairs for each base pair and distance interval
    for entry in data:
        base_pair, seq_i, seq_j, distance = entry
        for i, (lower, upper) in enumerate(distance_bins):
            if lower < distance <= upper:
                counts[base_pair][i] += 1
                total_counts[i] += 1
                total_pairs[base_pair] += 1
                break

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
    header = ["Base Pair", "Distance Interval", "Observed Frequency", "Reference Frequency", "Log-Ratio"]

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
    output_file_path = f"{base_name}_energy.txt"

    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')
        for row in tabular_data:
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")




# MAIN 

# STEP 1 : Process the PDB file and collect C3' data in tabular file
C3_tabular_file_path = extraction_C3_tabular_format(input_file_path)

# STEP 2 : Call the function to compute distances and write it in a tabular file
distance_tab_file = compute_distances(C3_tabular_file_path)

# STEP 3 : Calculate the estimated enregy for each base pairs and each distances intervals
calculate_energy(distance_tab_file)




