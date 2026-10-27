# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 15:55:34 2017

@author: TPflumm
"""
# Core Library modules
import math

# Third party modules
# THE PYTHON SHAPELY MODULE:
import shapely.geometry as shp_geom
# PythonOCC Libraries
from OCC.Core.BRepExtrema import BRepExtrema_DistShapeShape
from OCC.Core.GCPnts import GCPnts_QuasiUniformAbscissa
from OCC.Core.Geom2d import Geom2d_Line
from OCC.Core.Geom2dAdaptor import Geom2dAdaptor_Curve
from OCC.Core.Geom2dAPI import Geom2dAPI_InterCurveCurve
from OCC.Core.gp import (gp_Dir2d, gp_Lin2d, gp_Pnt2d, gp_Vec2d,)

# First party modules
from SONATA.cbm.mesh.mesh_utils import \
    remove_duplicates_from_list_preserving_order
from SONATA.cbm.topo.BSplineLst_utils import findPnt_on_curve


def map_node_on_curve(node, Curve2d, theta_11, distance=1e5, **kwargs):
    """
    maps a node onto a curve2d.

    Parameters
    ----------
    node : TYPE
        DESCRIPTION.
    Curve2d : TYPE
        DESCRIPTION.
    theta_11 : TYPE
        DESCRIPTION.
    distance : TYPE, optional
        DESCRIPTION. The default is 1e5.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """

    # ==================DIRECTION 1 ===============================
    Line1 = Geom2d_Line(gp_Lin2d(node.Pnt2d, gp_Dir2d(math.cos(math.radians(theta_11)), math.sin(math.radians(theta_11)))))
    # display.DisplayShape(Line1,color='BLACK')
    Intersection1 = Geom2dAPI_InterCurveCurve(Curve2d, Line1)
    dist = distance
    for i in range(1, Intersection1.NbPoints() + 1):
        if node.Pnt2d.Distance(Intersection1.Point(i)) < dist:
            dist = node.Pnt2d.Distance(Intersection1.Point(i))
            NewPnt1 = Intersection1.Point(i)

    # ===================DIRECTION 2 ===============================
    Line2 = Geom2d_Line(gp_Lin2d(node.Pnt2d, gp_Dir2d(-math.sin(math.radians(theta_11)), math.cos(math.radians(theta_11)))))
    # display.DisplayShape(Line2,color='WHITE')
    Intersection2 = Geom2dAPI_InterCurveCurve(Curve2d, Line2)
    dist = distance
    for i in range(1, Intersection2.NbPoints() + 1):
        if node.Pnt2d.Distance(Intersection2.Point(i)) < dist:
            dist = node.Pnt2d.Distance(Intersection2.Point(i))
            NewPnt2 = Intersection2.Point(i)

    # ==================CHOOSE NEW POINT===============================
    try:
        if node.Pnt2d.Distance(NewPnt1) > node.Pnt2d.Distance(NewPnt2):
            node.Pnt2d = NewPnt2
        else:
            node.Pnt2d = NewPnt1

    except NameError:
        try:
            node.Pnt2d = NewPnt2
        except NameError:
            node.Pnt2d = NewPnt1

    return None


def map_mesh_by_intersect_curve2d(mesh, curve2d, wire, global_minLen, **kwargs):
    """
    maps the mesh with and intersecting 2dcurve

    Parameters
    ----------
    mesh : list of cells
        DESCRIPTION.
    curve2d : curve2d
        DESCRIPTION.
    wire : OCC.topods.wire
        DESCRIPTION.
    global_minLen : TYPE
        DESCRIPTION.
    **kwargs : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    """
    # KWARGS:
    if kwargs.get("display") is not None:
        display = kwargs.get("display")
    else:
        display = None

    # # =========================INTERSECTING CELLS===============================
    icells = []
    for c in mesh:
        Distance = BRepExtrema_DistShapeShape(wire, c.wire)
        if Distance.Value() < 1e-5:
            icells.append(c)

    # Determine the interior nodes of Intersecting Cells
    # ===================INTERIOR NODES of Intersecting Cells===================
    Adaptor = Geom2dAdaptor_Curve(curve2d)
    NbPoints = 60
    discretization = GCPnts_QuasiUniformAbscissa(Adaptor, NbPoints)
    shpPointLst2 = []
    for j in range(1, NbPoints):
        para = discretization.Parameter(j)
        Pnt = gp_Pnt2d()
        curve2d.D0(para, Pnt)
        shpPointLst2.append((Pnt.X(), Pnt.Y()))

    shp_Poly = shp_geom.Polygon(shpPointLst2)

    inodes = []
    for c in icells:
        for n in c.nodes:
            shpPnt = shp_geom.Point(n.coordinates[0], n.coordinates[1])
            if shpPnt.within(shp_Poly):
                inodes.append(n)
                c.interior_nodes.append(n)

    inodes = sorted(set(inodes), key=lambda Node: (Node.id))
    #    for n in inodes:
    #        display.DisplayShape(n.Pnt2d,color='Orange')

    # =========================NONINTERSECTING INTERIOR CELLS===================
    niicells = []
    for c in mesh:
        inodes = []
        for n in c.nodes:
            shpPnt = shp_geom.Point(n.coordinates[0], n.coordinates[1])
            if shpPnt.within(shp_Poly):
                inodes.append(n)
        if len(inodes) == len(c.nodes):
            niicells.append(c)

    # for c in niicells:
    # display.DisplayShape(c.wire,color='RED')

    # =========================MAPPING==========================================
    mapped_nodes = []
    for c in icells:
        if len(c.nodes) == 4 and len(c.interior_nodes) == 1:
            node = c.interior_nodes[0]
            if node not in mapped_nodes:
                map_node_on_curve(node, curve2d, c.theta_11, display=display)
                mapped_nodes.append(node)

    for c in icells:
        if len(c.nodes) == 3 and len(c.interior_nodes) == 1:
            node = c.interior_nodes[0]
            if node not in mapped_nodes:
                map_node_on_curve(node, curve2d, c.theta_11, display=display)
                mapped_nodes.append(node)

    for c in icells:
        if len(c.nodes) == 4 and len(c.interior_nodes) == 2:
            node = c.interior_nodes[0]
            if node not in mapped_nodes:
                map_node_on_curve(node, curve2d, c.theta_11, display=display)
                mapped_nodes.append(node)

            node = c.interior_nodes[1]
            if node not in mapped_nodes:
                map_node_on_curve(node, curve2d, c.theta_11, display=display)
                mapped_nodes.append(node)

    rmcells = []
    for c in icells:
        if len(c.nodes) == 4 and len(c.interior_nodes) == 3:
            s = set(c.nodes) - set(c.interior_nodes)
            node = s.pop()
            if node not in mapped_nodes:
                map_node_on_curve(node, curve2d, c.theta_11, display=display)
                mapped_nodes.append(node)
            rmcells.append(c)

    for c in icells:
        if len(c.nodes) == 3 and len(c.interior_nodes) == 2:
            t = set(c.nodes) - set(mapped_nodes)
            if len(t) == 2:
                s = set(c.nodes) - set(c.interior_nodes)
                node = s.pop()

                if node not in mapped_nodes:
                    map_node_on_curve(node, curve2d, c.theta_11, display=display)
                    mapped_nodes.append(node)
                rmcells.append(c)

    # =========================POPPING CELLS====================================
    mesh = list(set(mesh) - set(niicells))
    mesh = list(set(mesh) - set(rmcells))

    rm2cells = []
    for c in list(set(mesh) - set(rmcells)):
        for n in c.nodes:
            shpPnt = shp_geom.Point(n.coordinates[0], n.coordinates[1])
            if shpPnt.within(shp_Poly):
                rm2cells.append(c)

    mesh = list(set(mesh) - set(rm2cells))

    # TODO: If two nodes are mapped to the same point, then
    # =========================SORTING MAPPED NODES=============================
    uLst = []
    mapped_nodes = list(set(mapped_nodes))
    for n in mapped_nodes:
        uLst.append(findPnt_on_curve(n.Pnt2d, curve2d))
    mapped_nodes = [x for y, x in sorted(zip(uLst, mapped_nodes))]

    # ====merge_mapped_nodes if to close:
    tmp = []
    for i, n1 in enumerate(mapped_nodes[0:], start=0):
        n2 = mapped_nodes[i - 1]
        v = gp_Vec2d(n1.Pnt2d, n2.Pnt2d)
        magnitude = v.Magnitude()

        if magnitude <= 0.001 * global_minLen:
            n1 = n2
        tmp.append(n1)

    mapped_nodes = remove_duplicates_from_list_preserving_order(tmp)

    return mesh, mapped_nodes
