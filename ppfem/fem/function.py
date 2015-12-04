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
    """
    Represents function on the given function space. The values can be initialized either via
      `set_dof_values_from_interpolation`
    or directly via
      `set_dof_values`.
    """
    def __init__(self, function_space):
        self.function_space = function_space
        self._dof_values = sp.zeros(self.function_space.number_of_dofs)

    def number_of_dofs(self):
        return self.function_space.number_of_dofs

    def function_dim(self):
        self.function_space.function_dim()

    def dof_values(self):
        return self._dof_values

    def set_dof_values_from_interpolation(self, function):
        """
        :param function: callable of kind f(x) where x is a vector of reference space dimension and
          and the return value if of dimension function_dim() of this function space
        """
        self.set_dof_values(self.function_space.interpolate_function(function))

    def set_dof_values(self, new_values):
        # TODO: add some checks
        self._dof_values[:] = new_values

    def localize(self, mesh_entity):
        elmt = self.function_space.get_element(mesh_entity)
        return LocalFEFunction(elmt,
                               self._dof_values[self.function_space.get_element_dof_index_array(elmt.index())])

    def get_element_dof_index_array(self, element_index):
        return self.function_space.get_element_dof_index_array(element_index)

    def get_mapping(self):
        return self.function_space.get_mapping()

    def get_mesh(self):
        return self.function_space.get_mesh()

    def get_subdomain(self):
        return self.function_space.get_subdomain()

    def __call__(self, mesh_entity, ref_point, der=0):
        if der == 0:
            return self.localize(mesh_entity).function_value(ref_point)
        elif der == 1:
            return self.localize(mesh_entity).function_gradient(ref_point)
        else:
            raise NotImplementedError("Derivatives higher than '1' are not implemented for FEFunction.")


class LocalFEFunction(object):
    """
    Class for element-level operations of FEFuntion. Objects should be instatiated by calling localize()
    on FEFunction objects.
    """
    def __init__(self, element, local_dof_values):
        self._element = element
        self._dof_values = local_dof_values

    def function_value(self, ref_point):
        """
        Computes the function value at a given point. If the underlying element is a mapped element,
        the point is expected to be in parameter space.
        :param ref_point: of type Point or a 1d array
        :return: a scalar of a numpy/scipy-array; depends on the underlying element.
        """
        return self._element.function_value(self._dof_values, ref_point)

    def function_gradient(self, ref_point):
        """
        Computes the function gradient at a given point. If the underlying element is a mapped element,
        the point is expected to be in parameter space.
        :param ref_point: of type Point or a 1d array
        :return: a numpy/scipy-array; depends on the underlying element.
        """
        return self._element.function_gradient(self._dof_values, ref_point)

    def __call__(self, ref_point, der=0):
        """
        This is simply a short-hand for function_value and function_gradient.
        See there for the details.
        :param der: "0" means value, "1" means gradient
        """
        if der == 0:
            return self.function_value(ref_point)
        elif der == 1:
            return self.function_gradient(ref_point)


class FunctionEvaluator(object):
    """
    This class acts as a wrapper for functions that are not an elements of a FE function space.
    Its only purpose is to provide an interface for function evaluations at points in parameter space such
    that the interface is similart to FEFunction. However, there is only a call operator implemented, no derivatives!
    Localization is provided following the general idea in PPFem.
    """
    def __init__(self, function, function_space):
        """
        Constructor
        :param function: a callable of kind f(x) where x is of type ppfem.geometry.point.Point
        :param function_space: an instance of FunctionSpace
        :return: an instance of FunctionEvaluator
        """
        self.function_space = function_space
        self.function = function

    def __call__(self, mesh_entity, ref_point):
        """
        Evaluates the wrapped function a the physical point that is obtained by mappping "ref_point" (an array
        or instance of Point) from the reference space to "phasical" space using the mapping between
        the given mesh_entity and the reference cell.
        :param mesh_entity: a mesh entity Line, Face, Cell
        :param ref_point: ana array or Point specifying a point in rerefence space
        :return: the return value of the wrapped function
        """
        return self.function_space.evaluate_function(self.function, mesh_entity, ref_point)

    def localize(self, mesh_entity):
        return LocalFunctionEvaluator(self.function, self.function_space.localize(mesh_entity))


class LocalFunctionEvaluator(object):
    """
    Class for element-level operations of FunctionEvaluator. Objects should be instatiated by calling localize()
    on FunctionEvaluator objects.
    """
    def __init__(self, function, local_function_space):
        self._function = function
        self._local_function_space = local_function_space

    def __call__(self, ref_point):
        """
        Evaluates the wrapped function a the physical point that is obtained by mappping "ref_point" (an array
        or instance of Point). Since a local function space is already associated to a mesh entity,
        ref_point is the only argument needed.
        :param ref_point: ana array or Point specifying a point in rerefence space
        :return: the return value of the wrapped function
        """
        return self._local_function_space.evaluate_function(self._function, ref_point)