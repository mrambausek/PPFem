from ppfem.fem.differential_equation import DifferentialEquation
import scipy as sp


class Laplace(DifferentialEquation):

    def __init__(self, dim):
        self._dim = dim

    def linear_form(self, element, quadrature, params):
        raise NotImplementedError("Implement me!")

    def bilinear_form(self, element, quadrature, params):
        raise NotImplementedError("Implement me!")
