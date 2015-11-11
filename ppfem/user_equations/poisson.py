from ppfem.user_equations.laplace import Laplace
import scipy as sp


class Poisson(Laplace):
    def __init__(self, rhs_function, dim):
        self._dim = dim
        self._rhs_function = rhs_function

    def compute_rhs(self, eval_data):
        return self._rhs_function(eval_data)
