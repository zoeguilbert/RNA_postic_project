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
    output_file_path = input_file_path.split('/')[-1].replace('.pdb', '_tabular.txt')
    with open(output_file_path, 'w') as output_file:
        output_file.write('\t'.join(header) + '\n')  # Write headers
        for row in tabular_data:                     # Write each parsed row
            output_file.write('\t'.join(map(str, row)) + '\n')

    print(f"Tabular file created: {output_file_path}")


# MAIN 
