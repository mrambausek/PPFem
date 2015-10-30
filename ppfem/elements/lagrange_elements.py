from ppfem.geometry.point import Point
from ppfem.elements.base import FiniteElement
from ppfem.elements.lagrange_basis import LagrangeBasis


class Line(FiniteElement):

    def __init__(self, degree, dimension=1):
        super(degree, dimension=dimension)

    def _setup_basis(self):
        support_points = self.get_support_points()
        self._basis_functions = [LagrangeBasis(support_points, i) for i in range(len(support_points))]

    def get_support_points(self):
        n = self._degree + 1
        return [Point(-1 + i * 2/n, index=i) for i in range(n)]








