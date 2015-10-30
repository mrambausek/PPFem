import numpy as np
import abc


class FiniteElement(abc.ABC):

    def __init__(self, degree, dimension):
        self._dimension = dimension
        self._degree = degree
        self._basis_functions = None
        self._setup_basis()

    @abc.abstractmethod
    def _setup_basis(self):
        raise Exception("Abstract method called!")

    def degree(self):
        return self._degree

    def dimension(self):
        return self._dimension

    def get_basis_function(self, index):
        return self._basis_functions[index]

    def get_basis_functions(self):
        return self._basis_functions

    def basis_function_value(self, index, point):
        return self._basis_functions[index].value(point)

    def basis_function_values(self, point):
        # first array axis corresponds to basis function!
        return np.array([f.value(point) for f in self._basis_functions])

    def basis_function_gradient(self, index, point):
        return self._basis_functions[index].gradient(point)

    def basis_function_gradients(self, point):
        # first array axis corresponds to basis function!
        return np.array([f.gradient(point) for f in self._basis_functions])

    def interpolate_value(self, dof_values, point):
        # first array axis corresponds to basis function!
        return np.tensordot(dof_values, self.basis_function_values(point), axes=([0], [0]))

    @abc.abstractmethod
    def get_support_points(self):
        raise Exception("Abstract method called!")
