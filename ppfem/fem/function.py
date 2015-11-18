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

import scipy as sp


class FEFunction(object):

    def __init__(self, function_space):
        self._function_space = function_space
        self._dof_values = sp.zeros(self._function_space.number_of_dofs)

    def number_of_dofs(self):
        return self._function_space.number_of_dofs

    def function_dim(self):
        self._function_space.function_dim()

    def dof_values(self):
        return self._dof_values

    def set_dof_values_from_interpolation(self, function):
        self.set_dof_values(self._function_space.interpolate_function(function))

    def set_dof_values(self, new_values):
        # TODO: add some checks
        self._dof_values = new_values

    def localize(self, mesh_entity):
        elmt = self._function_space.get_element(mesh_entity)
        return LocalFunction(elmt,
                             self._dof_values[self._function_space.get_element_dof_index_array(elmt.index())])

    def get_element_dof_index_array(self, element_index):
        return self._function_space.get_element_dof_index_array(element_index)


class LocalFunction(object):
    def __init__(self, element, local_dof_values):
        self._element = element
        self._dof_values = local_dof_values

    def function_value(self, ref_point):
        return self._element.function_value(self._dof_values, ref_point)

    def function_gradient(self, ref_point):
        return self._element.function_gradient(self._dof_values, ref_point)
