# PPFem: An educational finite element code
# Copyright (C) 2015  Matthias Rambausek
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import numpy as np
import abc


class ReferenceElement(abc.ABC):

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

    def basis_function_gradient(self, index, point, jacobian_inv=None):
        if jacobian_inv is None:
            return self._basis_functions[index].gradient(point)
        else:
            return np.dot(self._basis_functions[index].gradient(point), jacobian_inv)

    def basis_function_gradients(self, point, jacobian_inv=None):
        # first array axis corresponds to basis function!
        # return np.array([f.gradient(point) for f in self._basis_functions])
        return np.array([self.basis_function_gradient(i, point, jacobian_inv=jacobian_inv)
                         for i in range(self._n_bases)])

    @abc.abstractmethod
    def function_value(self, dof_values, point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def function_gradient(self, dof_values, point, jacobian_inv=None):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def get_support_points(self):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def interpolate_function(self, function, mapping=None):
        """
        :param function: callable accepting Point as argument return a numeric value or an array of numeric values
        :param mapping: optional, mapping of support points of the element from reference to physical space
        :return: an array of dof_values that in turn would reproduce the function values at the support points
        when given to fucntion_value(...)
        """
        raise Exception("Abstract method called!")

    @staticmethod
    @abc.abstractmethod
    def space_dim():
        raise Exception("Abstract method called!")

