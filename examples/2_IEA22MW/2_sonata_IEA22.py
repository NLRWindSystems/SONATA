import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval

# Path to yaml file
run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
job_str = 'IEA-22-280-RWT.yaml'
job_name = 'IEA22'
filename_str = run_dir + job_str

# ===== Define flags ===== #
flag_wt_ontology        = True # if true, use ontology definition of wind turbines for yaml files
flag_ref_axes_wt        = True # if true, rotate reference axes from wind definition to comply with SONATA (rotorcraft # definition)

# --- plotting flags ---
# Define mesh resolution, i.e. the number of points along the profile that is used for out-to-inboard meshing of a 2D blade cross section
mesh_resolution = 400
# For plots within blade_plot_sections
attribute_str           = 'MatID'  # default: 'MatID' (theta_3 - fiber orientation angle)
                                            # others:  'theta_3' - fiber orientation angle
                                            #          'stress.sigma11' (use sigma_ij to address specific component)
                                            #          'stressM.sigma11'
                                            #          'strain.epsilon11' (use epsilon_ij to address specific component)
                                            #          'strainM.epsilon11'

# 2D cross sectional plots (blade_plot_sections)
flag_plotTheta11        = False      # plane orientation angle
flag_recovery           = True
flag_plotDisplacement   = True     # Needs recovery flag to be activated - shows displacements from loadings in cross sectional plots
# 3D plots (blade_post_3dtopo)
flag_wf                 = True      # plot wire-frame
flag_lft                = True      # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create loft with grid refinement without too many radial_stations; can also export step file of lofted shape
flag_topo               = True      # plot mesh topology
c2_axis                 = False
flag_DeamDyn_def_transform = True               # transform from SONATA to BeamDyn coordinate system
flag_write_BeamDyn = True                       # write BeamDyn input files for follow-up OpenFAST analysis (requires flag_DeamDyn_def_transform = True)
flag_write_BeamDyn_unit_convert = ''  #'mm_to_m'     # applied only when exported to BeamDyn files
flag_estimate_structural_damping = True # Estimate structural damping from SONATA and BeamDyn
flag_output_zero_twist = False # Output sectional analysis so that BeamDyn sees everything at 0 twist?

# Shape of corners
choose_cutoff = 2    # 0 step, 2 round


# create flag dictionary
flags_dict = {"flag_wt_ontology": flag_wt_ontology, "flag_ref_axes_wt": flag_ref_axes_wt,
              "attribute_str": attribute_str,
              "flag_plotDisplacement": flag_plotDisplacement, "flag_plotTheta11": flag_plotTheta11,
              "flag_wf": flag_wf, "flag_lft": flag_lft, "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
              "flag_recovery": flag_recovery, "c2_axis": c2_axis}


# ===== User defined radial stations ===== #
# Define the radial stations for cross sectional analysis (only used for flag_wt_ontology = True -> otherwise, sections from yaml file are used!)

#radial_stations = [0.        , 0.02040816, 0.04081633, 0.06122449, 0.08163265, 0.10204082, 0.12244898, 0.14285714, 0.16326531, 0.18367347, 0.20408163, 0.2244898 , 0.24489796,
#               0.26530612, 0.28571429, 0.30612245, 0.32653061, 0.34693878, 0.36734694, 0.3877551 , 0.40816327, 0.42857143, 0.44897959, 0.46938776, 0.48979592, 0.51020408,
#               0.53061224, 0.55102041, 0.57142857, 0.59183673, 0.6122449 , 0.63265306, 0.65306122, 0.67346939, 0.69387755, 0.71428571, 0.73469388, 0.75510204, 0.7755102 ,
#               0.79591837, 0.81632653, 0.83673469, 0.85714286, 0.87755102, 0.89795918, 0.91836735, 0.93877551, 0.95918367, 1.]
radial_stations = np.array([0., 0.24489796, 0.51020408, 0.7755102, 1.])  # define radial stations as fraction of total blade length (from 0 to 1)

# ===== Execute SONATA Blade Component Object ===== #
# name          - job name of current task
# filename      - string combining the defined folder directory and the job name
# flags         - communicates flag dictionary (defined above)
# stations      - input of radial stations for cross sectional analysis
# stations_sine - input of radial stations for refinement (only and automatically applied when lofing flag flag_lft = True)
job = Blade(name=job_name, filename=filename_str, flags=flags_dict, cutoff_style = choose_cutoff, stations=radial_stations)  # initialize job with respective yaml input file

# ===== Build & mesh segments ===== #
job.blade_gen_section(topo_flag=True, mesh_flag = True)

# Define flags
flag_3d = False
flag_csv_export = True                         # export csv files with structural data
# Update flags dictionary
flags_dict['flag_csv_export'] = flag_csv_export
flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform
flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert
flags_dict['flag_output_zero_twist'] = flag_output_zero_twist
Loads_dict = {"Forces":[1.,1.,1.],"Moments":[1.,1.,1.]}

if not flag_estimate_structural_damping:
    # Set damping for BeamDyn input file
    delta = np.array([0.03, 0.03, 0.06787]) # logarithmic decrement, natural log of the ratio of the amplitudes of any two successive peaks. 3% flap and edge, 6% torsion
    zeta = 1. / np.sqrt(1.+(2.*np.pi / delta)**2.) # damping ratio,  dimensionless measure describing how oscillations in a system decay after a disturbance
    omega = np.array([0.508286, 0.694685, 4.084712])*2*np.pi # Frequency (rad/s), flap/edge/torsion
    mu1 = 2*zeta[0]/omega[0]
    mu2 = 2*zeta[1]/omega[1]
    mu3 = 2*zeta[2]/omega[2]
    mu = np.array([mu1, mu2, mu3, mu2, mu1, mu3])
else:
    mu = np.zeros(6)

beam_struct_eval(job_name, flags_dict, Loads_dict, radial_stations, job, run_dir,
                 job_str, mu)

if flag_estimate_structural_damping:
    # export stress-strain map for cross sectional analysis
    job.blade_exp_stress_strain_map(
        flag_output_zero_twist=flag_output_zero_twist)

# ===== PLOTS ===== #
# job.blade_plot_attributes()
# job.blade_plot_beam_props()

# saves figures in folder_str/figures if savepath is provided:
job.blade_plot_sections(attribute=attribute_str, plotTheta11=flag_plotTheta11,
                    plotDisplacement=flag_plotDisplacement, savepath=run_dir)
if flag_3d:
    job.blade_post_3dtopo(flag_wf=flags_dict['flag_wf'],
                          flag_lft=flags_dict['flag_lft'],
                          flag_topo=flags_dict['flag_topo'])
