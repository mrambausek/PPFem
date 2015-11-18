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


import abc


class LazyEval(object):
    def __init__(self, func, args):
        self._func = func
        self._args = args

    def __call__(self):
        return self._func(*self._args)


class PhysicalElement(abc.ABC):

    def __init__(self):
        pass
    # TODO: add more methods for extraction of local DoF indices for sub-entities
    # this might also make necessary storing related (additional) data in MeshEntity

    @abc.abstractmethod
    def number_of_global_dofs_per_vertex(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def number_of_global_non_vertex_dofs(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def number_of_global_dofs(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def topological_dimension(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def value_dimension(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def set_mesh_entity(self, mesh_entity):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def boundary_normal(self, local_boundary_index, boundary_ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def boundary_orientation(self, local_boundary_index):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def interpolate_function(self, mesh_entity):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def function_value(self, dof_values, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def function_gradient(self, dof_values, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def shape_function_values(self, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def shape_function_gradients(self, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def physical_coords(self, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def global_vertex_indices(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def index(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def domain_indicator(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def boundary_indicator(self):
        raise Exception('Abstract method called!')


class MappedElement(PhysicalElement):
    def __init__(self):
        PhysicalElement.__init__(self)
        self._mapping = None
        self._ref_element = None

    @abc.abstractmethod
    def set_mapping(self, mapping):
        raise Exception('Abstract method called!')

    def get_mapping(self):
        return self._mapping

    def get_reference_element(self):
        return self._ref_element
