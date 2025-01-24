"""
Author: Zoe GUILBERT
Date: 2025-01-24

Description:
This script loads tabular data from a file, categorizes nucleotide pairs, sorts distance intervals,
groups data by category, calculates mean log-ratios and observed frequencies, and visualizes the results.

Instructions:
1. Set the input file path to your tabular data file.
2. Run the script.

Input Table:
The input tabular file should contain the following columns:
- Base Pair: The base pair (e.g., "A-U").
- Distance Interval: The distance interval (e.g., "6-7").
- Observed Frequency: The observed frequency for the interaction.
- Log-Ratio: The log-ratio score for the interaction.

Output:
- Plots of mean log-ratios and observed frequencies for each nucleotide pair category.
"""

# Things to Change
"""
1. Set the input file path to your tabular data file.
"""
file_path = "/home/zozo/Documents/RNA_postic/02_results/combined_distances_energy.txt"



# Step 1: Load the tabular data
data = []

# Lire le fichier ligne par ligne
with open(file_path, 'r') as file:
    headers = file.readline().strip().split('\t')  # Extraire les en-têtes
    for line in file:
        row = line.strip().split('\t')
        row_dict = dict(zip(headers, row))
        row_dict['Observed_Frequency'] = float(row_dict['Observed_Frequency'])
        row_dict['Log-Ratio'] = float(row_dict['Log-Ratio'])
        row_dict['Distance_Interval'] = tuple(map(float, row_dict['Distance_Interval'].split('-')))  # Convertir intervalle en tuple
        data.append(row_dict)

# Step 2: Catégoriser les paires de bases
def nucleotide_type(base_pair):
    """Classifie le type de paire de bases."""
    big = {'A', 'G'}
    small = {'C', 'U'}
    b1, b2 = base_pair.split('-')
    if b1 in big and b2 in big:
        return "big-big"
    elif b1 in small and b2 in small:
        return "small-small"
    else:
        return "mixed"

# Ajouter une catégorie pour chaque ligne
for row in data:
    row['Category'] = nucleotide_type(row['Base_Pair'])

# Step 3: Trier les intervalles de distance
data.sort(key=lambda x: x['Distance_Interval'])  # Trier par distance réelle

# Step 4: Grouper les données par catégorie
from collections import defaultdict

grouped_log_ratios = defaultdict(list)
grouped_observed_frequencies = defaultdict(list)

for row in data:
    key = (row['Distance_Interval'], row['Category'])
    grouped_log_ratios[key].append(row['Log-Ratio'])
    grouped_observed_frequencies[key].append(row['Observed_Frequency'])

# Moyennes pour chaque groupe
mean_log_ratios = {key: sum(values) / len(values) for key, values in grouped_log_ratios.items()}
mean_observed_frequencies = {key: sum(values) / len(values) for key, values in grouped_observed_frequencies.items()}

# Préparer les données pour les plots
categories = ["big-big", "small-small", "mixed"]
x_labels = sorted(set(k[0] for k in grouped_log_ratios.keys()), key=lambda x: x[0])  # Trier les intervalles
x_positions = range(len(x_labels))

# Moyennes pour chaque catégorie
log_ratio_means = {cat: [mean_log_ratios.get((x, cat), 0) for x in x_labels] for cat in categories}
observed_freq_means = {cat: [mean_observed_frequencies.get((x, cat), 0) for x in x_labels] for cat in categories}

# Step 5: Visualisation
import matplotlib.pyplot as plt

# Plot des Log-Ratios
plt.figure(figsize=(10, 6))
for cat, means in log_ratio_means.items():
    plt.plot(x_positions, means, label=cat, marker='o')

plt.xticks(x_positions, [f"{int(x[0])}-{int(x[1])}" for x in x_labels], rotation=45)
plt.xlabel("Distance Interval")
plt.ylabel("Mean Log-Ratio")
plt.title("Log-Ratio Comparison of Nucleotide Pairs")
plt.legend()
plt.tight_layout()
plt.show()

# Barplot pour les fréquences observées
plt.figure(figsize=(10, 6))
width = 0.25  # Largeur des barres
for i, cat in enumerate(categories):
    means = observed_freq_means[cat]
    plt.bar(
        [x + i * width for x in x_positions],  # Décalage des barres
        means,
        width=width,
        label=cat
    )

plt.xticks(
    [x + width for x in x_positions],  # Centrer les étiquettes
    [f"{int(x[0])}-{int(x[1])}" for x in x_labels],
    rotation=45
)
plt.xlabel("Distance Interval")
plt.ylabel("Mean Observed Frequency")
plt.title("Observed Frequency Comparison of Nucleotide Pairs")
plt.legend()
plt.tight_layout()
plt.show()

# Main code
if __name__ == "__main__":
    # The script will run automatically when executed
    pass