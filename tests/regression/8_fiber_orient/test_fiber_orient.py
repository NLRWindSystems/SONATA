"""
Tests that verify fiber orientation behavior remains
consistent for future versions of SONATA.

Note: This test does not verify correctness,
just that the calculations of SONATA are stable.

The sign mapping flips for shear may be incorrect.
These are written as tests to verify SONATA is stable, not necessarily correct.

These seem to be a concern when they also require flipping the normal stresses/strains
(see mapping_signs_pm90)

"""

import os
import sys
import numpy as np
import pytest
import warnings

from SONATA.classBlade import Blade
from SONATA.utl.beam_struct_eval import beam_struct_eval

import matplotlib

matplotlib.use('Agg')

# Make the shared regression test utilities importable (located one directory up).
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))

import utils

def calculate_6x6_stress_maps(orientation, fiber_angle, split_quads=True):
    """
    With a fixed mesh, calculate the 6x6 matrices and stress/strain recovery maps.
    """


    # Path to yaml file
    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep
    job_str = './materials.yaml'

    # Encode the orientation, fiber angle and split_quads flag in the job name
    # so that the exported BeamDyn blade file is unique to this configuration.
    job_name = f'test-orient{orientation:}-fiber{fiber_angle:}-split{split_quads:}'
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
        theta_3 = fiber_angle

        for (x, cs) in job.sections:

            cs.cbm_custom_mesh(nodes, cells, MatID, split_quads=split_quads,
                            theta_11=theta_11, theta_3=theta_3)

    # ===== VABS / ANBAX combination/verification ===== #
    # Define flags
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

    beam_struct_eval(job_name, flags_dict, Loads_dict, radial_stations, job,
                    run_dir, job_str, mu)

    # ===== Output Recovery Matrices for Stress and Strain ===== #

    # The SONATA model currently has the unique stations.
    # Run stress maps for each of those stations.
    job.blade_exp_stress_strain_map(output_folder=os.path.join(run_dir, map_folder))
    map_fname = os.path.join(run_dir, map_folder, 'blade_station0000_stress_strain_map.npz')

    # BeamDyn blade file written by beam_struct_eval / write_beamdyn_prop as
    # os.path.join(folder_str, job_name + '_BeamDyn_Blade.dat'). Because job_name
    # now encodes the orientation, fiber angle and split_quads flag, this path is
    # unique to the current configuration.
    bd_fname = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')

    return bd_fname, map_fname

def compare_results(bd_fname1, map_fname1, bd_fname2, map_fname2,
                    expect_pass=True,
                    mass_rtol=1e-6, stiff_rtol=1e-6,
                    strain_rtol=1e-6, stress_rtol=1e-6,
                    mapping_signs=None):

    # Mapping signs flips the stress/strain map signs
    if mapping_signs is None:
        mapping_signs = np.ones(6)

    # ===== 6x6 mass and stiffness matrices ===== #
    mass1, stiff1 = utils.load_bd_blade(bd_fname1)
    mass2, stiff2 = utils.load_bd_blade(bd_fname2)

    mass1 = mass1[0]
    stiff1 = stiff1[0]
    mass2 = mass2[0]
    stiff2 = stiff2[0]

    assert mass1.shape == (6,6), 'Wrong shape of mass matrix.'
    assert stiff1.shape == (6,6), 'Wrong shape of stiffness matrix.'
    assert mass2.shape == (6,6), 'Wrong shape of mass matrix.'
    assert stiff2.shape == (6,6), 'Wrong shape of stiffness matrix.'

    mass_pass = np.abs(mass1 - mass2) < mass_rtol * np.abs(mass1).max()

    # always expect the mass to pass regardless of fiber orientations
    assert np.all(mass_pass), 'Mass matrices do not match.'


    stiff_pass = np.abs(stiff1 - stiff2) < stiff_rtol * np.abs(stiff1).max()

    if expect_pass:
        assert np.all(stiff_pass), 'Stiffness matrices do not match.'
    else:
        assert not np.all(stiff_pass), \
            'Stiffness matrices pass tolerances when they should not. Need to improve test.'


    # ===== Stress and strain recovery maps ===== #
    map1 = np.load(map_fname1)
    map2 = np.load(map_fname2)

    fc_to_strain_m_1 = map1['fc_to_strain_m']
    fc_to_stress_m_1 = map1['fc_to_stress_m']

    fc_to_strain_m_2 = map2['fc_to_strain_m'] * mapping_signs.reshape(1, -1, 1)
    fc_to_stress_m_2 = map2['fc_to_stress_m'] * mapping_signs.reshape(1, -1, 1)

    strain_pass = np.abs(fc_to_strain_m_1 - fc_to_strain_m_2) \
        < strain_rtol*np.abs(fc_to_strain_m_1).max()

    stress_pass = np.abs(fc_to_stress_m_1 - fc_to_stress_m_2) \
        < stress_rtol*np.abs(fc_to_stress_m_1).max()

    if expect_pass:
        assert np.all(strain_pass), 'Strain mapping matrices do not match.'
        assert np.all(stress_pass), 'Stress mapping matrices do not match.'
    else:

        assert not np.all(strain_pass), \
            'Strain mapping matrices pass tolerances when they should not. Need to improve test.'

        assert not np.all(stress_pass), \
            'Stress mapping matrices pass tolerances when they should not. Need to improve test.'

def test_orient0_fiber0():

    orientation = 0
    fiber_angle = 0
    split_quads = True # False does not appear to work anymore.

    bd_fname,map_fname = calculate_6x6_stress_maps(orientation,
                                                    fiber_angle,
                                                    split_quads=split_quads)

    #######################
    # Reference files

    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep

    job_name = f'./ref/test-orient{orientation:}-fiber{fiber_angle:}-split{True:}'

    bd_fname_ref = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation:}-fiber{fiber_angle:}-split{split_quads:}'

    map_fname_ref = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    #######################
    # Compare results

    compare_results(bd_fname_ref, map_fname_ref, bd_fname, map_fname,
                    expect_pass=True,
                    mass_rtol=1e-6, stiff_rtol=1e-6,
                    strain_rtol=1e-6, stress_rtol=1e-6)

    return

def test_orient0_fiber45():

    orientation = 0
    fiber_angle = 45
    split_quads = True # False does not appear to work anymore.

    bd_fname_p45,map_fname_p45 = calculate_6x6_stress_maps(orientation,
                                                        fiber_angle,
                                                        split_quads=split_quads)

    bd_fname_m45,map_fname_m45 = calculate_6x6_stress_maps(orientation,
                                                        -1*fiber_angle,
                                                        split_quads=split_quads)

    #######################
    # Reference files

    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep

    job_name = f'./ref/test-orient{orientation:}-fiber{fiber_angle:}-split{True:}'

    bd_fname_ref = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation:}-fiber{fiber_angle:}-split{split_quads:}'

    map_fname_ref = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    #######################
    # Compare results

    mass_rtol = 1e-6
    stiff_rtol = 1e-6
    strain_rtol = 1e-6
    stress_rtol = 1e-6

    # Match fiber angle response to reference.
    compare_results(bd_fname_ref, map_fname_ref, bd_fname_p45, map_fname_p45,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    # Test should be able to tell the difference between going the correct
    # and incorrect way for the fiber angle, at same tolerances
    compare_results(bd_fname_ref, map_fname_ref, bd_fname_m45, map_fname_m45,
                    expect_pass=False,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    return


def test_orient90_fiber0():

    orientation = 90
    fiber_angle = 0
    split_quads = True # False does not appear to work anymore.

    bd_fname_p90,map_fname_p90 = calculate_6x6_stress_maps(orientation,
                                                        fiber_angle,
                                                        split_quads=split_quads)

    bd_fname_m90,map_fname_m90 = calculate_6x6_stress_maps(-1*orientation,
                                                        fiber_angle,
                                                        split_quads=split_quads)

    # for stress/strain components these two configurations have the directions
    # flipped for 2 and 3
    # Right now, I am not convinced this is correct, but this should get
    # an exact match with ANBA.
    # Revist this later if tests stop passing.
    # this is what passes the test, but seems to be wrong.
    # mapping_signs_pm90 = np.array([1, -1, -1, -1, 1, 1])
    #
    # Should not need to flip any of the first three component signs.
    # Expect to need to flip the last two components, because the 1 direction is the same, but the
    # 2 and 3 are flipped so the shear only sees one sign flip
    # shear 23 sees two sign flips, so should be positive.
    mapping_signs_pm90 = np.array([1, 1, 1, 1, -1, -1])

    warnings.warn('These mapping sign flips are required for the test to pass, but seem wrong.'
        + ' Leaving this as a code stability test right now.')

    #######################
    # Reference files

    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep

    job_name = f'./ref/test-orient{orientation:}-fiber{fiber_angle:}-split{True:}'

    bd_fname_ref_90 = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation:}-fiber{fiber_angle:}-split{split_quads:}'

    map_fname_ref_90 = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')


    orientation_ref0 = 0

    job_name = f'./ref/test-orient{orientation_ref0:}-fiber{fiber_angle:}-split{True:}'

    bd_fname_ref_0 = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation_ref0:}-fiber{fiber_angle:}-split{split_quads:}'

    map_fname_ref_0 = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    #######################
    # Compare results

    mass_rtol = 1e-6
    stiff_rtol = 1e-6
    strain_rtol = 1e-6
    stress_rtol = 1e-6

    # Match orientation angle response to reference.
    compare_results(bd_fname_ref_90, map_fname_ref_90, bd_fname_p90, map_fname_p90,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    # match orientation +/-90 are the same
    compare_results(bd_fname_p90, map_fname_p90, bd_fname_m90, map_fname_m90,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol,
                    mapping_signs=mapping_signs_pm90)

    # Test should be able to tell the difference between orientation 0 and 90
    compare_results(bd_fname_ref_0, map_fname_ref_0, bd_fname_p90, map_fname_p90,
                    expect_pass=False,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    return


def test_orient90_fiber45():

    orientation = 90
    fiber_angle = 45
    split_quads = True # False does not appear to work anymore.

    bd_fname_p45,map_fname_p45 = calculate_6x6_stress_maps(orientation,
                                                        fiber_angle,
                                                        split_quads=split_quads)

    #######################
    # Reference files

    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep

    job_name = f'./ref/test-orient{orientation:}-fiber{fiber_angle:}-split{True:}'

    bd_fname_ref = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation:}-fiber{fiber_angle:}-split{split_quads:}'

    map_fname_ref = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    ###########
    # 0,0 reference
    fiber_angle_ref0 = 0
    orientation_ref0 = 0

    job_name = f'./ref/test-orient{orientation_ref0:}-fiber{fiber_angle_ref0:}-split{True:}'

    bd_fname_ref_0 = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation_ref0:}-fiber{fiber_angle_ref0:}-split{split_quads:}'

    map_fname_ref_0 = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    #######################
    # Compare results

    mass_rtol = 1e-6
    stiff_rtol = 1e-6
    strain_rtol = 1e-6
    stress_rtol = 1e-6

    # Match fiber angle response to reference.
    compare_results(bd_fname_ref, map_fname_ref, bd_fname_p45, map_fname_p45,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    # Test should be able to tell the difference between going the correct
    # and incorrect way for the fiber angle, at same tolerances
    compare_results(bd_fname_ref_0, map_fname_ref_0, bd_fname_p45, map_fname_p45,
                    expect_pass=False,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    return


def test_orient40_fiber20():

    orientation_p = 40
    fiber_angle_p = 20

    orientation_m = -140
    fiber_angle_m = -20

    split_quads = True # False does not appear to work anymore.

    bd_fname_pp,map_fname_pp = calculate_6x6_stress_maps(orientation_p,
                                                       fiber_angle_p,
                                                       split_quads=split_quads)

    bd_fname_pm,map_fname_pm = calculate_6x6_stress_maps(orientation_p,
                                                       fiber_angle_m,
                                                       split_quads=split_quads)

    bd_fname_mp,map_fname_mp = calculate_6x6_stress_maps(orientation_m,
                                                       fiber_angle_p,
                                                       split_quads=split_quads)

    bd_fname_mm,map_fname_mm = calculate_6x6_stress_maps(orientation_m,
                                                       fiber_angle_m,
                                                       split_quads=split_quads)

    # for stress/strain components these two configurations have the directions
    # flipped for 2 and 3
    # Right now, I am not convinced this is correct, but this should get
    # an exact match with ANBA.
    # Revist this later if tests stop passing.
    #
    # Should not need to flip any of the first three component signs.
    # Expect to need to flip the last two components, because the 1 direction is the same, but the
    # 2 and 3 are flipped so the shear only sees one sign flip
    # shear 23 sees two sign flips, so should be positive.
    mapping_signs = np.array([1, 1, 1, 1, -1, -1])

    warnings.warn('These mapping sign flips are required for the test to pass, but seem wrong.'
        + ' Leaving this as a code stability test right now.')

    #######################
    # Reference files


    #################
    # +/+ reference

    run_dir = os.path.dirname( os.path.realpath(__file__) ) + os.sep

    job_name = f'./ref/test-orient{orientation_p:}-fiber{fiber_angle_p:}-split{True:}'

    bd_fname_ref_pp = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation_p:}-fiber{fiber_angle_p:}-split{split_quads:}'

    map_fname_ref_pp = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    #################
    # +/- reference

    job_name = f'./ref/test-orient{orientation_p:}-fiber{fiber_angle_m:}-split{True:}'

    bd_fname_ref_pm = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation_p:}-fiber{fiber_angle_m:}-split{split_quads:}'

    map_fname_ref_pm = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')


    #################
    # 0,0 reference

    orientation_ref0 = 0
    fiber_angle_ref0 = 0

    job_name = f'./ref/test-orient{orientation_ref0:}-fiber{fiber_angle_ref0:}-split{True:}'

    bd_fname_ref_0 = os.path.join(run_dir, job_name + '_BeamDyn_Blade.dat')
    map_folder = f'stress-maps-orient{orientation_ref0:}-fiber{fiber_angle_ref0:}-split{split_quads:}'

    map_fname_ref_0 = os.path.join(run_dir, 'ref', map_folder, 'blade_station0000_stress_strain_map.npz')

    #######################
    # Compare results

    mass_rtol = 1e-6
    stiff_rtol = 1e-6
    strain_rtol = 1e-6
    stress_rtol = 1e-6

    # Match orientation angle response to reference.
    compare_results(bd_fname_ref_pp, map_fname_ref_pp, bd_fname_pp, map_fname_pp,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    # Flipping both angles should get back to the same results
    compare_results(bd_fname_pp, map_fname_pp, bd_fname_mm, map_fname_mm,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol,
                    mapping_signs=mapping_signs)


    # Reference check for flipping only one angle
    compare_results(bd_fname_ref_pm, map_fname_ref_pm, bd_fname_pm, map_fname_pm,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    # Results should match regardless of which angle is flipped
    compare_results(bd_fname_pm, map_fname_pm, bd_fname_mp, map_fname_mp,
                    expect_pass=True,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol,
                    mapping_signs=mapping_signs)

    # Results should be different between positive and negative fiber angles
    compare_results(bd_fname_pp, map_fname_pp, bd_fname_pm, map_fname_pm,
                    expect_pass=False,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    # Test should be able to tell the difference between orientation 0 and 90
    compare_results(bd_fname_ref_0, map_fname_ref_0, bd_fname_pp, map_fname_pp,
                    expect_pass=False,
                    mass_rtol=mass_rtol, stiff_rtol=stiff_rtol,
                    strain_rtol=strain_rtol, stress_rtol=stress_rtol)

    return



if __name__ == "__main__":
    pytest.main(["-s", "test_fiber_orient.py"])
