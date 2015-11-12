from ppfem.fem.partial_differential_equation import PDE
import scipy as sp


class Laplace(PDE):

    def __init__(self, test_function_space, trial_function_space, quadrature, dim=1):
        PDE.__init__(self, test_function_space, trial_function_space, quadrature)
        self._dim = dim

    def local_linear_form(self, params):
        raise NotImplementedError("Implement me!")

    def local_bilinear_form(self, params):
        raise NotImplementedError("Implement me!")
