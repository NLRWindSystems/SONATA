import matplotlib.pyplot as plt
import numpy as np
import os
import yaml

# Load data from the file
run_dir = os.path.dirname(os.path.abspath(__file__))
file_path_damp = os.path.join(run_dir, "damp_material_contributions.txt")
file_path_energy = os.path.join(run_dir, "energy_material_contributions.txt")
file_modal = os.path.join(run_dir, "modal_data.yaml")

damp_contributions = np.loadtxt(file_path_damp)
energy_contributions = np.loadtxt(file_path_energy)
mat_names = ["UD GFRP", "BX GFRP", "TX GFRP", "Adhesive", "Foam", "Gelcoat"]
Mod_names = ["1F", "1E", "2F", "2E", "3F", "3E", "4F", "1T"]
with open(file_modal, 'r') as f:
    modal_data = yaml.safe_load(f)

Freqs = [modal_data['freq_Hz'][i] for i in range(len(Mod_names))]

for i in range(len(Mod_names)):
    # Create a pie chart for each mode
    labels = mat_names[:]
    # Lump slices smaller than 3% into a single slice
    threshold = 3  # percentage threshold

    damp_data = damp_contributions[i, :]
    damp_small_slices = damp_data < (threshold / 100) * np.sum(damp_data)
    if np.any(damp_small_slices):
        damp_data = np.append(damp_data[~damp_small_slices], np.sum(damp_data[damp_small_slices]))
        damp_labels = [labels[j] for j in range(len(labels)) if not damp_small_slices[j]] + ['Other']

    energy_data = energy_contributions[i, :]
    energy_small_slices = energy_data < (threshold / 100) * np.sum(energy_data)
    if np.any(energy_small_slices):
        energy_data = np.append(energy_data[~energy_small_slices], np.sum(energy_data[energy_small_slices]))
        energy_labels = [labels[j] for j in range(len(labels)) if not energy_small_slices[j]] + ['Other']

    plt.figure(figsize=(12, 6))
    colors = plt.cm.tab20.colors  # Use a colormap with enough distinct colors
    color_map = {label: colors[i % len(colors)] for i, label in enumerate(mat_names + ['Other'])}
    slice_colors = [color_map[label] for label in labels]
    plt.title(f'Mode {Mod_names[i]} - Freq: {Freqs[i]} Hz ', weight='bold')
    plt.axis('off')


    plt.subplot(1, 2, 1)
    plt.pie(energy_data, labels=energy_labels, autopct='%1.1f%%', textprops={'weight': 'bold'}, colors=[color_map[label] for label in energy_labels])
    plt.title('Strain Energy', weight='bold')

    plt.subplot(1, 2, 2)
    plt.pie(damp_data, labels=damp_labels, autopct='%1.1f%%', textprops={'weight': 'bold'}, colors=[color_map[label] for label in damp_labels])
    plt.title('Damping', weight='bold')
    # plt.pie(data, labels=labels, autopct='%1.1f%%', textprops={'weight': 'bold'}, colors=slice_colors)
    # plt.pie(contributions[i, :], labels=mat_names, autopct='%1.1f%%', startangle=140, textprops={'weight': 'bold'})
    plt.savefig(os.path.join(run_dir, f'material_contributions_mode_{Mod_names[i]}.png'), dpi=300)
    plt.savefig(os.path.join(run_dir, f'material_contributions_mode_{Mod_names[i]}.pdf'), dpi=300)
    plt.tight_layout()
    plt.close()
