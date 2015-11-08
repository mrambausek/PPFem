from ppfem.fem.material import Material
import scipy as sp


class Laplace(Material):
    def __init__(self, dim):
        self._dim = dim

    def compute_rhs(self, physical_point, function_value, function_gradient):
        return sp.zeros((self._dim,1))

    def compute_lhs(self, physical_point, function_value, function_gradient):
        return sp.eye(self._dim)