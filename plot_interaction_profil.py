"""
Author: Zoe GUILBERT
Date: 2025-01-24

Description:
This script reads a tabular file containing interaction score for each nucleotide pairs.
It also plots the interaction profiles for each base pair, and saves the plots in a specified folder.

Instructions:
1. Set the input file path to your tabular file containing interaction data.
2. Set the output folder path where you want to save the plots.
3. Run the script.

Input Table:
The input tabular file should contain interaction data with the following columns:
- Base Pair: The base pair (e.g., "A-U").
- Distance Interval: The distance interval (e.g., "6-7").
- Log-Ratio: The log-ratio score for the interaction.

Output:
- Plots of interaction profiles for each base pair saved in the specified output folder.
"""


import matplotlib.pyplot as plt
import os

# Define the input file and output folder where you want to save the plots
input_file = '/home/zozo/Documents/RNA_postic/02_results_rna_prediction/combined_distances_energy.txt'
output_folder = "/home/zozo/Documents/RNA_postic/02_results_rna_prediction/interaction_profiles"

# Create the folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Save the plot in the specified folder

# Function to plot the interaction profiles
def plot_interaction_profiles(tabular_file_path, output_folder):
    # Read the tabular file and parse the data
    with open(tabular_file_path, 'r') as file:
        lines = file.readlines()

    # Skip the header
    data = [line.strip().split('\t') for line in lines[1:]]

    # Extract unique base pairs
    base_pairs = set(entry[0] for entry in data)

    # Create a plot for each base pair
    for base_pair in base_pairs:
        # Filter the data for the current base pair
        pair_data = [entry for entry in data if entry[0] == base_pair]

        # Extract the distance intervals and log-ratios
        distances = [float(interval.split('-')[0]) + 0.5 for interval in [entry[1] for entry in pair_data]]
        log_ratios = [float(entry[4]) for entry in pair_data]

        # Plot the interaction profile
        plt.figure(figsize=(10, 6))
        plt.plot(distances, log_ratios, marker='o', linestyle='-', label=base_pair)
        plt.xlabel('Distance (Å)')
        plt.ylabel('Log-Ratio (Score)')
        plt.title(f'Interaction Profile for {base_pair}')
        plt.grid(True)
        plt.legend()
        plt.savefig(f"{output_folder}/interaction_profile_{base_pair.replace('-', '_')}.png")
        plt.close()

    print("Interaction profiles plotted and saved in", output_folder)


#Main code
if __name__ == "__main__":
    # Call the function to plot interaction profiles
    plot_interaction_profiles(input_file, output_folder)
