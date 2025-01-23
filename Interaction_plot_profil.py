import matplotlib.pyplot as plt
import os

# Define the folder where you want to save the plots
output_folder = "/home/zozo/Documents/RNA_postic/RNA_postic_project/interaction_profiles"

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

    print("Interaction profiles plotted and saved.")

# Call the function to plot interaction profiles
plot_interaction_profiles('/home/zozo/Documents/RNA_postic/RNA_postic_project/1gid_C3_tabular_distances_energy_tabular.txt', output_folder)
