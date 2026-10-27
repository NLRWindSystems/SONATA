"""
Use OpenFAST linearization to:
    - Calculate mode shapes
    - Plot mode shapes
    - Calculate modal internal forces and moments
    - Save a YAML file of the internal forces/moments

BeamDyn should not have any damping in the simulation


"""

import os
import numpy as np
from openfast_toolbox.io.fast_linearization_file import FASTLinearizationFile

import yaml
import matplotlib
import matplotlib.pyplot as plt

scriptDir = os.path.dirname(__file__)

###############################################################################
######## Inputs (including argparse)                                   ########
###############################################################################

#####################
# Folders

run_dir = os.path.dirname(os.path.realpath(__file__))
map_folder = os.path.join(os.path.dirname(run_dir), "stress-map")
show_plots = 1

if not show_plots:
    matplotlib.use('Agg')


#####################

# lin_fname = './openfast/Free_Free_Blade_Linear.1.BD1.lin'
lin_fname = os.path.join(run_dir,
                         'openfast/IEA-22-280-RWT.1.lin')

bd_sum_fname = os.path.join(run_dir,
                         'openfast/IEA-22-280-RWT.BD1.sum.yaml')

output_fname = os.path.join(run_dir,
                         'modal_data.yaml')

output_fname_energy = os.path.join(run_dir,
                         'modal_strain_energy_data.yaml')



N_modes = 10

###############################################################################
######## openfast toolbox example                                      ########
###############################################################################

print('To do free free analysis, you need to augment the modes with the tower base motion and proceed from there.')

# --- Open lin File
# Files:
#    Ideal_Beam_Fixed_Free_Linear.1.BD1.lin
#    Ideal_Beam_Fixed_Free_Linear.1.ED.lin
#    Ideal_Beam_Fixed_Free_Linear.1.lin
linFile = os.path.join(scriptDir, lin_fname)

print('Line File:')
print(linFile)

# linFile = os.path.join(scriptDir, 'Fixed_Free_Blade_Linear.1.lin')
lin = FASTLinearizationFile(linFile)
# print('Keys available:',lin.keys())

# --- Perform eigenvalue analysis
#fd, zeta, Q, f0 = eigA(lin['A'])
fd, zeta, mode_shapes, f0 = lin.eva()
print('Nat. freq. [Hz], Damping ratio [%]')
print(np.column_stack((np.around(f0,4),np.around(zeta*100,4))))

# --- Using dataframes instead of numpy arrays for convenient indexing of variables
dfs = lin.toDataFrame()
# print('Dataframe available:',dfs.keys())
A = dfs['A']
# print(A.columns)


###############################################################################
######## Manual Eigenvalue Analysis to Verify 6 RBM                    ########
###############################################################################

v = np.linalg.eig(A)[0]

v = v[np.argsort(np.abs(v))]

omega_0 = np.abs(v)              # natural cylic frequency [rad/s]
zeta_np    = - np.real(v)/omega_0   # damping ratio
freq_0  = omega_0/(2*np.pi)      # natural frequency [Hz]

print('Manual Calc [Freq [Hz], Freq [Hz], zeta [%], zeta [%]]')

print(np.hstack((freq_0.reshape(-1, 2), zeta_np.reshape(-1, 2))))

freq_short = np.round(freq_0[:12+2*N_modes:2], 10)

freq_str = ''.join('{:.2f}, '.format(val) for val in freq_short)[:-2]

print('Short Freq list [Hz]: ' + freq_str)


###############################################################################
######## Eliminate Rigid Body Modes                                    ########
###############################################################################

# Doubled modes from manual eigenanalysis, 6 rigid body modes expected
# assert np.max(np.abs(freq_0[:12])) / freq_0[12] < 1e-2

# first_freq = freq_0[12]

# mask_no_rbm = f0 >= 0.5 * first_freq

# fd = fd[mask_no_rbm]
# zeta = zeta[mask_no_rbm]
# mode_shapes = mode_shapes[:, mask_no_rbm]
# f0 = f0[mask_no_rbm]

###############################################################################
######## Load BD Summary File                                          ########
###############################################################################

with open(bd_sum_fname) as f:
    data = list(yaml.load_all(f, Loader=yaml.loader.SafeLoader))[-1]

# Node initial position and rotation
node_x0 = np.array(data['Init_Nodes_E1'])[:,2]
quads_x0 = np.array(data['Init_QP_E1'])[:, 2]

###############################################################################
######## Create Masks for DOFs to Plot Mode Shapes                     ########
###############################################################################

dof_masks = {'Flap' : [False] * len(A.columns),
             'Edge' : [False] * len(A.columns),
             'Axial' : [False] * len(A.columns),
             'EdgeR' : [False] * len(A.columns),
             'FlapR' : [False] * len(A.columns),
             'Torsion' : [False] * len(A.columns)}

# # Templates for just the BD linearization file
# dof_name_template = {'Flap' : '^N.*FLAP-D',
#                      'Edge' : '^N.*EDGE-D',
#                      'Axial' : '^N.*transZ',
#                      'EdgeR' : '^N.*EDGE-R',
#                      'FlapR' : '^N.*FLAP-R',
#                      'Torsion' : '^N.*TORS-R'}

# Templates for using the full openfast
dof_name_template = {'Flap' : '^BD_B1N.*FLAP-D',
                     'Edge' : '^BD_B1N.*EDGE-D',
                     'Axial' : '^BD_B1N.*transZ',
                     'EdgeR' : '^BD_B1N.*FLAP-R',
                     'FlapR' : '^BD_B1N.*EDGE-R',
                     'Torsion' : '^BD_B1N.*TORS-R'}

for key in dof_masks.keys():
    dof_masks[key] = A.keys().str.contains(dof_name_template[key])

# A.keys()[dof_masks['Flap']]

###############################################################################
######## Plot Mode Shapes                                              ########
###############################################################################

for dof in dof_masks.keys():

    # Plot Displacements
    for i in range(N_modes):

        angle = np.angle(np.hstack(([0], mode_shapes[dof_masks[dof],  i])))

        plt.plot(node_x0,
                 np.hstack(([0], np.abs(mode_shapes[dof_masks[dof],  i]))) \
                     * np.sign(angle),
                 label='Mode {}'.format(i+1))

    plt.title(dof + ' Abs*Sign(Angle)')
    plt.legend()
    # plt.show()

    # # Plot Phase Angles
    # for i in range(N_modes):
    #     angle = np.angle(np.hstack(([0], mode_shapes[dof_masks[dof],  i])))
    #     plt.plot(node_x0, angle,
    #              label='Mode {}'.format(i+1))
    # plt.title(dof + ' Phase Angle')
    # plt.legend()
    # plt.show()

###############################################################################
######## Recover Internal Forces/Moments                               ########
###############################################################################

C = dfs['C']

for ind,key in enumerate(A.columns):
    assert C.columns[ind] == key, 'Ordering of A and C do not match.'


output_data = {}

norm_ref = 1e-12

for key in ['FxL', 'FyL', 'FzL', 'MxL', 'MyL', 'MzL']:

    mask = C.index.str.contains(key)

    internal = C.loc[mask, :].values @ mode_shapes[:, :N_modes]

    norm_ref = np.maximum(norm_ref, np.linalg.norm(internal))

    # Check that everything is nearly completely imaginary
    assert np.linalg.norm(np.real(internal)) / norm_ref < 1e-8, \
        'Internal forces are not predominately imaginary as assumed here.' \
        + ' Make sure there is no damping in the simulation.'

    output_data[key] = np.imag(internal)

###############################################################################
######## Plot Internal Forces / Moments                                ########
###############################################################################

for key in output_data.keys():

    # Plot Displacements
    for i in range(N_modes):

        plt.plot(quads_x0, output_data[key][:, i],
                 label='Mode {}'.format(i+1))

    plt.title(key)
    plt.legend()
    # plt.show()

###############################################################################
######## Add Displacement Info to Outputs                              ########
###############################################################################

output_data['Nodes_Z'] = node_x0
output_data['Quads_Z'] = quads_x0

output_data['tip_flap'] = np.abs(mode_shapes[dof_masks['Flap'], :][-1, :N_modes]
                                 ).tolist()

output_data['tip_edge'] = np.abs(mode_shapes[dof_masks['Edge'], :][-1, :N_modes]
                                 ).tolist()

output_data['tip_twist'] = np.abs(mode_shapes[dof_masks['Torsion'], :][-1, :N_modes]
                                  ).tolist()

output_data['freq_Hz'] = freq_short.tolist()

###############################################################################
######## Save outputs                                                  ########
###############################################################################

# Convert all numpy arrays to lists before saving
for key in output_data.keys():
    if isinstance(output_data[key],np.ndarray):
        output_data[key] = output_data[key].tolist()

    if isinstance(output_data[key], list):
        if isinstance(output_data[key][0], list):
            # Nested list
            for i in range(len(output_data[key])):
                output_data[key][i] = [float(val) for val in output_data[key][i]]
        else:
            # Single list only
            output_data[key] = [float(val) for val in output_data[key]]
    else:
        # Scalar
        output_data[key] = float(output_data[key])

# Write out the yaml file of the modal information
with open(output_fname, 'w') as file:
    yaml.dump(output_data, file)


###############################################################################
######## Load Some OTTR Mode Data to Run                               ########
###############################################################################

# N_modes = 1

# csv_file = '/Users/jporter/repos/rust-bd/ottr/output/bar-subcomponent/' \
#     + 'rot_000_rad_s/fc_star_eigen_07.csv'

# csv_data = np.loadtxt(csv_file, delimiter=',', skiprows=2)

# output_data['FzL'] = csv_data[:, 2:3]
# output_data['FyL'] = csv_data[:, 3:4]
# output_data['FxL'] = -1*csv_data[:, 4:5]
# output_data['MzL'] = csv_data[:, 5:6]
# output_data['MyL'] = csv_data[:, 6:7]
# output_data['MxL'] = -1*csv_data[:, 7:8]

###############################################################################
######## Stress and Strain Energy Recovery                             ########
###############################################################################

map_data = [None] * len(output_data['FzL'])

for i in range(len(output_data['FzL'])):
    map_fname = os.path.join(map_folder, 'blade_station{:04d}_stress_strain_map.npz'.format(i))

    map_data[i] = np.load(map_fname)

# quadrature locations
quad_pos = np.array(output_data['Quads_Z'])

# weights for trapezoid integration between stations
quad_length_weights = np.zeros(len(quad_pos))
quad_length_weights[0]    = 0.5*(quad_pos[1] - quad_pos[0])
quad_length_weights[-1]   = 0.5*(quad_pos[-1] - quad_pos[-2])
quad_length_weights[1:-1] = 0.5*(quad_pos[2:] - quad_pos[0:-2])

energy_dict = {}
for station_ind, map_curr in enumerate(map_data):
    for mat_ind, mat_name in enumerate(map_curr['material_names']):
        energy_dict[mat_name] = np.zeros((6, N_modes))


# # omit the None material and the gelcoat
# for material in map_data['material_names'][1:-1]:
#     energy_dict[str(material)] = np.full((6, N_modes), 0.)

# Loop over modes
for mode_ind in range(N_modes):

    force_moments = np.vstack((np.array(output_data['FzL'])[:, mode_ind],
                              np.array(output_data['FyL'])[:, mode_ind],
                              -np.array(output_data['FxL'])[:, mode_ind],
                              np.array(output_data['MzL'])[:, mode_ind],
                              np.array(output_data['MyL'])[:, mode_ind],
                              -np.array(output_data['MxL'])[:, mode_ind]))

    for station_ind, map_curr in enumerate(map_data):
        # After einsum:
        # First index is stress/strain component in order [11, 22, 33, 23, 13, 12]
        # second index is element number
        stresses = np.einsum('ijk,j->ik', map_curr['fc_to_stress_m'],
                            force_moments[:, station_ind])

        strains = np.einsum('ijk,j->ik', map_curr['fc_to_strain_m'],
                            force_moments[:, station_ind])

        # energy densities at each element and each station
        energy_densities = 0.5*stresses * strains

        # energy scaled by volume for each element (area and length quadrature
        # weight)
        energy = energy_densities * map_curr['elem_areas'].reshape(1, -1) \
                    * quad_length_weights[station_ind]

        # loop over materials and add up energies for each direction
        for mat_ind, mat_name in enumerate(map_curr['material_names']):

            elem_mask = map_curr['elem_materials'] == mat_ind

            if elem_mask.sum() > 0:

                # first index for each material is direction with order
                # [11, 22, 33, 23, 13, 12]


                # sum energy over elements and then over stations
                energy_dict[mat_name][:, mode_ind] \
                    += energy[:, elem_mask].sum(axis=1)

###############################################################################
######## For a given mode, print energy fractions                      ########
###############################################################################

mode_ind = 0

mode_energy = np.zeros(N_modes)


for mode_ind in range(N_modes):
    print('Mode {}:'.format(mode_ind + 1))
    for key in energy_dict.keys():

        mode_energy[mode_ind] += energy_dict[key][:, mode_ind].sum()

    for key in energy_dict.keys():
        print('{} Energy Fractions: {}'.format(key,
                       energy_dict[key][:, mode_ind]/mode_energy[mode_ind]))

        energy_dict[key][:, mode_ind] = energy_dict[key][:, mode_ind] \
                                            / mode_energy[mode_ind]


    # Replace the energy dictionary with the fraction for each mode


###############################################################################
######## Save strain energy outputs                                    ########
###############################################################################

# Convert all numpy arrays to lists before saving
for key in energy_dict.keys():
    if isinstance(energy_dict[key],np.ndarray):
        energy_dict[key] = energy_dict[key].tolist()

    if isinstance(energy_dict[key], list):
        if isinstance(energy_dict[key][0], list):
            # Nested list
            for i in range(len(energy_dict[key])):
                energy_dict[key][i] = [float(val) for val in energy_dict[key][i]]
        else:
            # Single list only
            energy_dict[key] = [float(val) for val in energy_dict[key]]
    else:
        # Scalar
        energy_dict[key] = float(energy_dict[key])

# Write out the yaml file of the modal information
with open(output_fname_energy, 'w') as file:
    yaml.dump(energy_dict, file)


# Order is [11, 22, 33, 23, 13, 12]
loss_factors = {'glass_uniax' : np.array([1.49e-3, 7.16e-3, 7.16e-3,
                                       6.62e-3, 6.62e-3, 6.62e-3]),
                'carbon_uniax' : np.array([4.69e-4, 5.61e-3, 5.61e-3,
                                       5.61e-3, 5.61e-3, 5.61e-3]),
                'glass_biax' : np.array([7.04e-3, 7.04e-3, 7.16e-3,
                                      6.62e-3, 6.62e-3, 3.01e-3]),
                'glass_triax' : np.array([4.05E-03, 7.00E-03, 7.16E-03,
                                       6.62E-03, 6.62E-03, 4.82E-03]),
                'Adhesive' : np.array([1.3e-2]*6),
                'medium_density_foam' : np.array([1.3e-2]*6), # Same of adhesive
                'gelcoat' : np.array([1.3e-2]*6), # same of gelcoat
                }

for key in energy_dict.keys():
    energy_dict[key] = np.array(energy_dict[key])

###############################################################################
######## Load 1D + 2D Loss Factors and Damping                         ########
###############################################################################

loss_factors_vec = np.hstack((loss_factors['glass_uniax'],
                              loss_factors['carbon_uniax'],
                              loss_factors['glass_biax'],
                              loss_factors['glass_triax'],
                              loss_factors['Adhesive'],
                              loss_factors['medium_density_foam'],
                              loss_factors['gelcoat']))


energy_fracs_mat = np.hstack((energy_dict['glass_uniax'].T,
                              energy_dict['carbon_uniax'].T,
                              energy_dict['glass_biax'].T,
                              energy_dict['glass_triax'].T,
                              energy_dict['Adhesive'].T,
                              energy_dict['medium_density_foam'].T,
                              energy_dict['gelcoat'].T))



sonata_zeta = 0.5 * (energy_fracs_mat @ loss_factors_vec)

print("Critical damping ratios (%):")
print(sonata_zeta*100.)
critical_damping_ratios_file = os.path.join(run_dir, 'critical_damping_ratios.txt')
np.savetxt(critical_damping_ratios_file, sonata_zeta*100.)


damp_contributions = np.zeros((N_modes, len(loss_factors)))
energy_contributions = np.zeros((N_modes, len(loss_factors)))
for i in range(N_modes):
    for j, material in enumerate(loss_factors.keys()):
        energy_contributions[i,j] = sum(energy_dict[material][:, i])
        damp_contributions[i,j] = sum(energy_dict[material][:, i] * loss_factors[material])

# Save contributions and labels to a text file
damp_contributions_file = os.path.join(run_dir, 'damp_material_contributions.txt')
energy_contributions_file = os.path.join(run_dir, 'energy_material_contributions.txt')

np.savetxt(energy_contributions_file, energy_contributions)
np.savetxt(damp_contributions_file, damp_contributions)

# Save loss factors as a table
loss_factors_table_file = os.path.join(run_dir, 'loss_factors_table.txt')

with open(loss_factors_table_file, 'w') as file:
    file.write("Material\tLoss Factors (11, 22, 33, 23, 13, 12)\n")
    for material, factors in loss_factors.items():
        factors_str = ', '.join(f'{factor:.6e}' for factor in factors)
        file.write(f"{material}\t{factors_str}\n")
