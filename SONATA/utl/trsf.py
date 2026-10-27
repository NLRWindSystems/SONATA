#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 09:48:14 2019

@author: gu32kij
"""
# Third party modules
import numpy as np
from OCC.Core.gp import (gp_Ax1, gp_Ax2, gp_Ax3, gp_Dir,
                         gp_Pnt, gp_Trsf,)


def trsf_blfr_to_cbm(Ax2_blfr, Ax2_befr):
    """
    generates the transform from the SONATA reference frame to the local 2d cbm
    frame

    Parameters
    --------
    Ax2_blfr : gp_Ax2
        Opencascade blade reference frame,
        Describes a right-handed coordinate system in 3D space.
    Ax2_befr : gp_Ax2
        Opencascade locale beam reference frame,
        Opencascade: Describes a right-handed coordinate system in 3D space.

    Returns
    --------
    Trsf : gp_Trsf
        Opencascade: non-persistent transformation in 3D space

    """
    Ax2_cbm = gp_Ax2(Ax2_befr.Location(), Ax2_befr.XDirection(), Ax2_befr.YDirection())
    # print('Ax2_cbm.Direction:', Ax2_cbm.XDirection().Coord(), Ax2_cbm.YDirection().Coord())
    trsf1 = gp_Trsf()
    trsf2 = gp_Trsf()

    trsf1.SetTransformation(gp_Ax3(Ax2_blfr), gp_Ax3(Ax2_befr))
    trsf2.SetTransformation(gp_Ax3(Ax2_befr), gp_Ax3(Ax2_cbm))

    Trsf = gp_Trsf()
    Trsf.Multiply(trsf2)
    Trsf.Multiply(trsf1)

    return Trsf


def trsf_cbm_to_blfr(Ax2_blfr, Ax2_befr):
    """
    generates the transform from the cbm local 2d frame to the blade frame

    Parameters
    --------
    Ax2_blfr : gp_Ax2
        Opencascade blade reference frame,
        Describes a right-handed coordinate system in 3D space.
    Ax2_befr : gp_Ax2
        Opencascade locale beam reference frame,
        Opencascade: Describes a right-handed coordinate system in 3D space.

    Returns
    --------
    Trsf : gp_Trsf
        Opencascade: non-persistent transformation in 3D space

    """
    Ax2_cbm = gp_Ax2(Ax2_befr.Location(), Ax2_befr.XDirection(), Ax2_befr.YDirection())
    # print('Ax2_cbm.Direction:', Ax2_cbm.XDirection().Coord(), Ax2_cbm.YDirection().Coord())

    trsf1 = gp_Trsf()
    trsf2 = gp_Trsf()

    trsf1.SetTransformation(gp_Ax3(Ax2_cbm), gp_Ax3(Ax2_befr))
    trsf2.SetTransformation(gp_Ax3(Ax2_befr), gp_Ax3(Ax2_blfr))

    Trsf = gp_Trsf()
    Trsf.Multiply(trsf2)
    Trsf.Multiply(trsf1)

    return Trsf


def trsf_af_to_blfr(loc, soy, chord, twist, deformation=None):
    """
    Defines the transformation in 3D space to the blade reference frame location
    and pitch-axis information, scales it with chord information and rotates
    it with twist information

    Parameters
    ----------
    loc : array
        [x,y,z] position in blade reference coordinates
    soy : float
        dim. pitch axis location
    chord : float
        chordlength
    twist : float
        twist angle about x in radians

    Returns
    ---------
    Trsf : OCC.gp_Trsf
        non-persistent transformation in 3D space.

    Todo
    ---------
    implement to transfer to the deformed stated, check for rotational definition!
        [xnew,ynew,znew,phi,theta,psi] deformation vector is defined as...


    """
    trsf_rot1 = gp_Trsf()
    trsf_rot2 = gp_Trsf()
    trsf_rot3 = gp_Trsf()
    trsf_trans1 = gp_Trsf()
    trsf_trans2 = gp_Trsf()
    trsf_scale = gp_Trsf()

    pa_loc = soy / chord

    trsf_rot1.SetRotation(gp_Ax1(gp_Pnt(pa_loc, 0, 0), gp_Dir(0, 0, 1)), -twist)
    trsf_trans1.SetTranslation(gp_Pnt(pa_loc, 0, 0), gp_Pnt(0, 0, 0))
    trsf_scale.SetScale(gp_Pnt(0, 0, 0), chord)
    trsf_rot2.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 0, 1)), -np.pi / 2)
    trsf_rot3.SetRotation(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(0, 1, 0)), -np.pi / 2)
    trsf_trans2.SetTranslation(gp_Pnt(0, 0, 0), gp_Pnt(loc[0], loc[1], loc[2]))

    if deformation:
        trsf_trans2 = gp_Trsf()
        trsf_deform_rot1 = gp_Trsf()
        trsf_deform_rot2 = gp_Trsf()
        trsf_deform_rot3 = gp_Trsf()

        trsf_trans2.SetTranslation(gp_Pnt(0, 0, 0), gp_Pnt(deformation[0], deformation[1], deformation[2]))
        trsf_deform_rot1.SetRotation(gp_Ax1(gp_Pnt(loc[0], loc[1], loc[2]), gp_Dir(1, 0, 0)), deformation[3])
        trsf_deform_rot2.SetRotation(gp_Ax1(gp_Pnt(loc[0], loc[1], loc[2]), gp_Dir(0, 1, 0)), deformation[4])
        trsf_deform_rot3.SetRotation(gp_Ax1(gp_Pnt(loc[0], loc[1], loc[2]), gp_Dir(0, 0, 1)), deformation[5])

    Trsf = gp_Trsf()
    Trsf.Multiply(trsf_trans2)
    Trsf.Multiply(trsf_rot3)
    Trsf.Multiply(trsf_rot2)
    Trsf.Multiply(trsf_scale)
    Trsf.Multiply(trsf_trans1)
    Trsf.Multiply(trsf_rot1)
    if deformation:
        Trsf.Multipy(trsf_deform_rot1)
        Trsf.Multipy(trsf_deform_rot2)
        Trsf.Multipy(trsf_deform_rot3)
    return Trsf


if __name__ == "__main__":
    Ax2_1 = gp_Ax2()
    Ax2_2 = gp_Ax2()

    trfs1 = trsf_cbm_to_blfr(Ax2_1, Ax2_2)
    trfs2 = trsf_blfr_to_cbm(Ax2_1, Ax2_2)
