import numpy as np
import abc


class FiniteElement(abc.ABC):
    # todo: information if a dof is internal or not

    def __init__(self, degree, dimension):
        self._dimension = dimension
        self._degree = degree
        self._basis_functions = None
        self._n_dofs = None
        self._n_internal_dofs = None
        self._n_bases = None
        self._setup_basis()

    @abc.abstractmethod
    def _setup_basis(self):
        """
        setup basis functions and set the number of dofs and basis functions
        """
        raise Exception("Abstract method called!")

    def degree(self):
        return self._degree

    def n_of_dofs(self):
        return self._n_dofs

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
        # return np.array([f.value(point) for f in self._basis_functions])
        return np.array([self.basis_function_value(i, point) for i in range(self._n_bases)])

    def basis_function_gradient(self, index, point):
        return self._basis_functions[index].gradient(point)

    def basis_function_gradients(self, point):
        # first array axis corresponds to basis function!
        # return np.array([f.gradient(point) for f in self._basis_functions])
        return np.array([self.basis_function_gradient(i, point) for i in range(self._n_bases)])

    def function_value(self, dof_values, point):
        # first array axis corresponds to basis function!
        if self._dimension == 1:
            # return np.float64(np.einsum('i,i', dof_values, self.basis_function_values(point)))
            return np.dot(dof_values, self.basis_function_values(point))
        else:
            return np.einsum('ijk,ijk->jk', dof_values, self.basis_function_values(point))

    def function_gradient(self, dof_values, point):
        # first array axis corresponds to basis function!
        if self._dimension == 1:
            # return np.float64(np.einsum('i,i', dof_values, self.basis_function_values(point)))
            return np.dot(dof_values, self.basis_function_gradients(point))
        elif self.space_dim() > 1:
            return np.einsum('ijk,ijkl->jkl', dof_values, self.basis_function_gradients(point))
        elif self.space_dim() == 1:
            return np.einsum('ijk,ijk->jk', dof_values, self.basis_function_gradients(point))

    @abc.abstractmethod
    def get_support_points(self):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def interpolate_function(self, function, mapping=None):
        """
        :param function: callable accepting Point as argument return a numeric value or an array of
        :param mapping: optional, mapping of support points of the element from reference to physical space
        :return: an array of dof_values that in turn would reproduce the function values at the support points
        when given to fucntion_value(...)
        """
        raise Exception("Abstract method called!")

    @staticmethod
    @abc.abstractmethod
    def space_dim():
        raise Exception("Abstract method called!")

