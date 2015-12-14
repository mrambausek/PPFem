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

from ppfem.geometry.point import Point
from ppfem.elements.base import ReferenceElement
from ppfem.elements.lagrange_basis import LagrangeBasis
import numpy as np


class LagrangeElement(ReferenceElement):
    def __init__(self, degree, dimension=1):
        ReferenceElement.__init__(self, degree, dimension)

    def interpolate_function(self, function, mapping=None):
        """
        This implementation shows the characteristic property of Lagrange Elements!
        :param function: a callable f(p) where p is a Point (or a coordinate array of size space_dim()) and the result
        is of dimension dimension()
        :param mapping: a Mapping instance to compute the "physical" coordinates of a point in reference space
        """
        if mapping is not None:
            points = mapping.map_points(self.get_support_points())
        else:
            points = self.get_support_points()

        return np.array([function(p.coords()) for p in points])

    def function_value(self, dof_values, point):
        # first array axis corresponds to basis function!
        if self._dimension == 1:
            return np.dot(self.basis_function_values(point).reshape(1, self._n_bases), dof_values)
        else:
            return np.einsum('ijk,ijk->j',
                             dof_values.reshape(self._n_bases, self._dimension, 1),
                             self.basis_function_values(point))

    def function_gradient(self, dof_values, point, jacobian_inv=None):
        # first array axis corresponds to basis function!
        if self._dimension == 1:
            return np.dot(self.basis_function_gradients(point, jacobian_inv=jacobian_inv).reshape(dof_values.shape).T,
                          dof_values)
        elif self.space_dim() > 1:
            # FIXME: shape problems
            print('1', dof_values.shape)
            print('2', self.basis_function_gradients(point, jacobian_inv=jacobian_inv).shape)
            raise NotImplementedError("Fix the shapes!")
            return np.einsum('ijk,ijkl->jkl',
                             dof_values,
                             self.basis_function_gradients(point, jacobian_inv=jacobian_inv))
        elif self.space_dim() == 1:
            return np.einsum('ijk,ijk->jk',
                             dof_values.reshape(self._n_bases, self._dimension, 1),
                             self.basis_function_gradients(point, jacobian_inv=jacobian_inv))


class LagrangeLine(LagrangeElement):
    def __init__(self, degree, dimension=1):
        LagrangeElement.__init__(self, degree, dimension=dimension)

    def _setup_basis(self):
        support_points = self.get_support_points()
        self._n_bases = len(support_points)
        self._n_dofs = self._n_bases * self._dimension
        self._n_internal_dofs = self._n_dofs - 2
        self._basis_functions = [LagrangeBasis(support_points, i, dimension=self._dimension)
                                 for i in range(len(support_points))]

    def get_support_points(self):
        n = self._degree + 1
        return [Point(-1), Point(1)] + [Point(-1 + i * 2 / (n - 1), index=i) for i in range(1, n - 1)]

    @staticmethod
    def space_dim():
        return 1
