from ppfem.user_materials.laplace import Laplace
import scipy as sp


class Poisson(Laplace):
    def __init__(self, rhs_function, dim):
        self._dim = dim
        self._rhs_function = rhs_function

    def compute_rhs(self, physical_point, function_value, function_gradient):
        return self._rhs_function(physical_point, function_value, function_gradient)
