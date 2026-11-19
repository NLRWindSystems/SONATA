from types import SimpleNamespace
from unittest.mock import Mock, patch

from SONATA.cbm.mesh import mesh_utils


def test_grab_nodes_on_bspline_list_skips_projection_domain_error():
    node = Mock(Pnt2d=object(), parameters=[0, 0, 0])
    nodes = [node]
    projection = SimpleNamespace(
        NbPoints=lambda: 1,
        Distance=lambda _: 0,
        LowerDistanceParameter=lambda: 0.25,
    )

    with patch.object(
        mesh_utils,
        "Geom2dAPI_ProjectPointOnCurve",
        side_effect=[ValueError("OpenCASCADE domain error"), projection],
    ):
        result = mesh_utils.grab_nodes_on_BSplineLst(nodes, [object(), object()])

    assert result == nodes
    assert node.parameters == [0, 1, 0.25]
