# PPFem: A educational finite element code
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

    def __call__(self, *args, **kwargs):
        return self._func(*self._args)


class PhysicalElement(abc.ABC):

    def __init__(self):
        pass

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
    def shape_function_gradients(self, dof_values, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def physical_coords(self, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def jxw(self, qp_data):
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
