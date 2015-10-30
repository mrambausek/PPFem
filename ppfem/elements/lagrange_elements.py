from ppfem.geometry.point import Point
from ppfem.elements.base import FiniteElement
from ppfem.elements.lagrange_basis import LagrangeBasis
import numpy as np


class LagrangeElement(FiniteElement):

    def __init__(self, degree, dimension=1):
        FiniteElement.__init__(self, degree, dimension)

    def interpolate_function(self, function, mapping=None):
        """
        This implementation shows the characteristic property of Lagrange Elements!
        """
        if mapping is not None:
            points = mapping.map_points(self.get_support_points())
        else:
            points = self.get_support_points()

        return np.array([function(p) for p in points])


class Line(LagrangeElement):

    def __init__(self, degree, dimension=1):
        LagrangeElement.__init__(self, degree, dimension=dimension)

    def _setup_basis(self):
        support_points = self.get_support_points()
        self._n_bases = len(support_points)
        self._n_dofs = self._n_bases * self._dimension
        self._basis_functions = [LagrangeBasis(support_points, i, dimension=self._dimension)
                                    for i in range(len(support_points))]

    def get_support_points(self):
        n = self._degree + 1
        return [Point(-1 + i * 2/(n-1), index=i) for i in range(n)]

    @staticmethod
    def space_dim():
        return 1
