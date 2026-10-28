"""
Tests that verify fiber orientation behavior remains
consistent for future versions of SONATA.

Note: This test does not verify correctness,
just that the calculations of SONATA are stable.

"""

import os
import numpy as np
from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval
# from SONATA.utl_openfast.utl_run_fast import utl_run_fast
from openfast_toolbox.io import FASTInputFile
import yaml

import mesh_routines

import matplotlib

matplotlib.use('Agg')

def calculate_6x6_stress_maps(orientation, fiber_angle, split_quads=True):
    """
    With a fixed mesh, calculate the 6x6 matrices and stress/strain recovery maps.
    """


    # Path to yaml file
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_str = './materials.yaml.yaml'

    job_name = 'test'
    filename_str = run_dir + job_str

    # Options for recovering stress/strain fields for modal response
    recover_stress = False # running stress maps, not individual calcs.
    custom_mesh = True

    #####################

    map_folder = f'stress-maps-orient{orientation:}-fiber{fiber_angle:}-split{split_quads:}'

    # ===== Define flags ===== #
    # if true, use ontology definition of wind turbines for yaml files
    flag_wt_ontology        = True

    # if true, rotate reference axes from wind definition to comply with SONATA
    #  (rotorcraft # definition)
    flag_ref_axes_wt        = True

    # --- plotting flags ---
    # Define mesh resolution, i.e. the number of points along the profile that is
    # used for out-to-inboard meshing of a 2D blade cross section
    mesh_resolution = 400

    # For plots within blade_plot_sections
    attribute_str           = 'theta_11'  # default: 'MatID' (theta_3 - fiber orientation angle)
                                    # others:  'theta_3' - fiber orientation angle
                                    #          'stress.sigma11' (use sigma_ij to address specific component)
                                    #          'stressM.sigma11'
                                    #          'strain.epsilon11' (use epsilon_ij to address specific component)
                                    #          'strainM.epsilon11'

    # 2D cross sectional plots (blade_plot_sections)
    flag_plotTheta11        = False      # plane orientation angle
    flag_recovery           = recover_stress

    # Needs recovery flag to be activated - shows displacements from loadings in
    #  cross sectional plots
    flag_plotDisplacement   = True

    # 3D plots (blade_post_3dtopo)
    flag_wf                 = True # plot wire-frame

    # plot lofted shape of blade surface (flag_wf=True obligatory); Note: create
    # loft with grid refinement without too many radial_stations; can also export
    # step file of lofted shape
    flag_lft                = True

    # plot mesh topology
    flag_topo               = True
    c2_axis                 = False


    # create flag dictionary
    flags_dict = {"flag_wt_ontology": flag_wt_ontology,
                "flag_ref_axes_wt": flag_ref_axes_wt,
                "attribute_str": attribute_str,
                "flag_plotDisplacement": flag_plotDisplacement,
                "flag_plotTheta11": flag_plotTheta11,
                "flag_wf": flag_wf, "flag_lft": flag_lft,
                "flag_topo": flag_topo, "mesh_resolution": mesh_resolution,
                "flag_recovery": flag_recovery, "c2_axis": c2_axis}


    # ===== User defined radial stations ===== #
    # Define the radial stations for cross sectional analysis (only used for
    # flag_wt_ontology = True -> otherwise, sections from yaml file are used!)
    radial_stations = [0.0]


    # ===== Execute SONATA Blade Component Object ===== #
    # name          - job name of current task
    # filename      - string combining the defined folder directory and the job name
    # flags         - communicates flag dictionary (defined above)
    # stations      - input of radial stations for cross sectional analysis
    # stations_sine - input of radial stations for refinement (only and
    #                 automatically applied when lofing flag flag_lft = True)

    # initialize job with respective yaml input file
    job = Blade(name=job_name, filename=filename_str, flags=flags_dict,
                stations=radial_stations)

    if custom_mesh:

        # 3x3 grid of nodes, x and y values of [-0.05, 0.0, 0.05]
        #   6 --- 7 --- 8      y = 0.05
        #   |     |     |
        #   3 --- 4 --- 5      y = 0.0
        #   |     |     |
        #   0 --- 1 --- 2      y = -0.05
        nodes = np.array([
            [-0.05, -0.05],  # 0
            [ 0.00, -0.05],  # 1
            [ 0.05, -0.05],  # 2
            [-0.05,  0.00],  # 3
            [ 0.00,  0.00],  # 4
            [ 0.05,  0.00],  # 5
            [-0.05,  0.05],  # 6
            [ 0.00,  0.05],  # 7
            [ 0.05,  0.05],  # 8
        ])

        # 4 quad elements, each listed clockwise starting from its bottom-left node:
        # bottom-left -> top-left -> top-right -> bottom-right
        cells = np.array([
            [0, 3, 4, 1],  # lower-left
            [1, 4, 5, 2],  # lower-right
            [3, 6, 7, 4],  # upper-left
            [4, 7, 8, 5],  # upper-right
        ])

        MatID = np.array([1, 1, 1, 1])
        theta_11 = np.ones(4) * orientation
        theta_3 = np.ones(4) * fiber_angle

        for (x, cs) in job.sections:

            cs.cbm_custom_mesh(nodes, cells, MatID, split_quads=split_quads,
                            theta_11=theta_11, theta_3=theta_3)

    # ===== VABS / ANBAX combination/verification ===== #
    # Define flags
    flag_3d = False
    flag_run_vabs = False
    flag_run_anbax = True

    # needs flag_run_vabs & flag_run_anbax set to True to be valid!
    flag_verify_vabs_anbax = False

    # transform from SONATA to BeamDyn coordinate system
    flag_DeamDyn_def_transform = True

    # plots and saves figures of structural characteristics of respective analysis
    # for the range of defined radial stations
    flag_plot_vabs_struct_characteristics = False

    # export csv files with structural data
    flag_csv_export = False

    # write BeamDyn input files for follow-up OpenFAST analysis
    # (requires flag_DeamDyn_def_transform = True)
    flag_write_BeamDyn = True

    # applied only when exported to BeamDyn files
    flag_write_BeamDyn_unit_convert = ''  #'mm_to_m'
    # Update flags dictionary
    flags_dict['flag_run_vabs'] = flag_run_vabs
    flags_dict['flag_run_anbax'] = flag_run_anbax
    flags_dict['flag_verify_vabs_anbax'] = flag_verify_vabs_anbax
    flags_dict['flag_DeamDyn_def_transform'] = flag_DeamDyn_def_transform

    flags_dict['flag_plot_vabs_struct_characteristics'] \
        = flag_plot_vabs_struct_characteristics

    flags_dict['flag_csv_export'] = flag_csv_export
    flags_dict['flag_write_BeamDyn'] = flag_write_BeamDyn
    flags_dict['flag_write_BeamDyn_unit_convert'] = flag_write_BeamDyn_unit_convert

    Loads_dict = {}

    # Set damping for BeamDyn input file
    mu = np.zeros((6,))

    beam_struct_eval(flags_dict, Loads_dict, radial_stations, job, run_dir,
                    job_str, mu)

    # ===== Output Recovery Matrices for Stress and Strain ===== #

    # The SONATA model currently has the unique stations.
    # Run stress maps for each of those stations.
    job.blade_exp_stress_strain_map(output_folder=map_folder)
    map_fname = os.path.join(map_folder, 'blade_station0000_stress_strain_map.npz')

    """
    Here, need to determine the BeamDyn blade filename that has the information for the
    6x6 matrices. This also means, above job probably needs to be modified so that
    the beamdyn blade file name is defined to include the orientation and fiber angle.
    """


    return bd_fname, map_fname

def compare_results(bd_fname1, map_fname1, bd_fname2, map_fname2):

    pass

def test_orient0_fiber0_triangles():

    pass