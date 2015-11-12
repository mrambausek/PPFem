from ppfem.fem.differential_equation import PDE
import scipy as sp


class Laplace(PDE):

    def __init__(self, element, quadrature, dim=1):
        PDE.__init__(self, element, quadrature)
        self._dim = dim

    def linear_form(self, mesh_entity, dof_values, params):
        raise NotImplementedError("Implement me!")

    def bilinear_form(self, mesh_entity, dof_values, params):
        raise NotImplementedError("Implement me!")
