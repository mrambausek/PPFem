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
from ppfem.fem.physical_element import MappedElement


class FunctionSpace(object):
    """
    This class represents FEM-style function spaces.
    They rely on a mesh, an element describing the local shape functions. Some additional input like a mapping object
    or a subdomain to which this function space is restrict can also be provided.
    FunctionSpace provide assembly information for elements via
      `get_element_dof_index_array`
    and allow to interpolate functions via
      `interpolate_function`
    where the result are the values of the values for "degrees of freedom" for an interpolant on this function space.
    """

    def __init__(self, element, mesh=None, subdomain=None, mapping=None):
        self._element = element
        self._mesh = mesh
        self._subdomain = subdomain
        self._element_dof_map = {}
        self._vertex_dof_map = {}
        self.number_of_dofs = None

        self.storage_ready = False
        if mesh is not None:
            self._setup_storage()

        if mapping is not None:
            self._element.set_mapping(mapping)

    def _setup_storage(self):
        for v in self._mesh.vertices():
            self._vertex_dof_map[v.global_index()] = []

        self.number_of_dofs = self._generate_assembly_data()
        self.storage_ready = True

    def set_mesh(self, mesh, sub_domain=None):
        self._mesh = mesh
        self._subdomain = sub_domain
        self._setup_storage()

    def get_mapping(self):
        if isinstance(self._element, MappedElement):
            return self._element.get_mapping()
        else:
            return None

    def get_mesh(self):
        return self._mesh

    def get_subdomain(self):
        return self._subdomain

    def mesh_entity_iterator(self, topological_dim=None):
        if self._subdomain is None:
            return self._mesh.get_mesh_entities(topological_dim=topological_dim)
        return filter(lambda e: e.domain_indicator == self._subdomain, self._mesh.get_mesh_entities())

    def _generate_element_vertex_dofs(self, element, first_dof_index, visited_vertices):
        vertex_dofs_per_element = element.number_of_global_dofs_per_vertex()

        dofs = []
        for v in element.global_vertex_indices():
            if v not in visited_vertices:
                next_first_dof_index = first_dof_index + vertex_dofs_per_element
                new_dofs = range(first_dof_index, next_first_dof_index)
                first_dof_index = next_first_dof_index
                self._vertex_dof_map[v] += new_dofs
            dofs += self._vertex_dof_map[v]

        return dofs, first_dof_index

    def _generate_element_non_vertex_dofs(self, element, first_dof_index):
        non_vertex_dofs_per_element = element.number_of_global_non_vertex_dofs()
        next_first_dof_index = first_dof_index+non_vertex_dofs_per_element
        dofs = range(first_dof_index, first_dof_index+non_vertex_dofs_per_element)
        return dofs, next_first_dof_index

    def _generate_assembly_data(self):
        """Generates and tabulates the global dofs for each element."""
        visited_vertices = []
        first_dof_index = 0

        for e in self._mesh.get_mesh_entities():
            if self._subdomain is not None and e.domain_indicator != self._subdomain:
                continue

            element = self.get_element(e)

            dofs = []
            new_dofs, first_dof_index = self._generate_element_vertex_dofs(
                element, first_dof_index, visited_vertices
            )
            dofs += new_dofs
            visited_vertices += e.global_vertex_indices()

            new_dofs, first_dof_index = self._generate_element_non_vertex_dofs(element, first_dof_index)
            dofs += new_dofs

            self._element_dof_map[element.index()] = dofs
        return first_dof_index

    def get_element_dof_index_array(self, element_index):
        return sp.array(self._element_dof_map[element_index], dtype=sp.int64)

    def get_element(self, mesh_entity):
        if self._subdomain is not None and mesh_entity.domain_indicator != self._subdomain:
            raise Exception("Mesh entity in wrong subdomain!")
        self._element.set_mesh_entity(mesh_entity)
        return self._element

    def get_number_of_global_element_dofs(self, mesh_entity):
        return self.get_element(mesh_entity).number_of_global_dofs()

    def localize(self, mesh_entity):
        return LocalFunctionSpace(self.get_element(mesh_entity))

    def function_dim(self):
        self._element.value_dimension()

    def function_gradient_shape(self):
        return self._element.value_dimension(), self._element.topological_dimension()

    def function_space_dim(self):
        return self.number_of_dofs

    def interpolate_function(self, function):
        """
        Computes the values for degrees of freedom in this function space to represent an interpolation of
        the given function.
        :param function: callable of kind f(x) where x is a vector of reference space dimension and
          and the return value if of dimension function_dim() of this function space
        :return: an array of values of degrees of freedom, which may be used for construction an object of type
        FEFunction
        """
        global_dofs = sp.zeros(self.number_of_dofs)
        for e in self.mesh_entity_iterator():
            global_dofs[self._element_dof_map[e.index]] = self.localize(e).interpolate_function(function)
        return global_dofs

    def evaluate_function(self, function, mesh_entity, ref_point):
        """
        The function value for 'ref_point' in 'mesh_entity'.
        :param function: callable of kind f(x) where x is a vector of physical space dimension and
          and the return value if of dimension function_dim() of this function space
        """
        return self.localize(mesh_entity).evaluate_function(function, ref_point)


class LocalFunctionSpace(object):
    def __init__(self, element):
        self._element = element

    def interpolate_function(self, function):
        return self._element.interpolate_function(function)

    def shape_function_values(self, ref_point):
        """
        Computes the values of the shape functions at a given point.
        For mapped elements, this point is expected to be in parameter space.
        :param ref_point: of type Point or a 1d-array
        :return: an array with the first axis correspondig to the local index of the shape function.
        The remaining shape is the shape of the function value. This depends on the underlying element.
        """
        return self._element.shape_function_values(ref_point)

    def shape_function_gradients(self, ref_point):
        """
        Computes the values of the shape functions at a given point.
        For mapped elements, this point is expected to be in parameter space.
        :param ref_point: of type Point or a 1d-array
        :return: an array with the first axis correspondig to the local index of the shape function.
        The remaining shape is the shape of the gradient of a shape function. This depends on the underlying element.
        (numpy.array or scipy.array respectively)
        """
        return self._element.shape_function_gradients(ref_point)

    def physical_coords(self, ref_point):
        return self._element.physical_coords(ref_point)

    def number_of_shape_functions(self):
        """
        The number of shape functions (from the underlying element). Shape functions may be
        vector-valued, so function_dim() may be also useful for you.
        :return: number of shape function (integer)
        """
        return self._element.number_of_global_dofs()

    def function_dim(self):
        """
        :return: value dimenstion of shape functions (integer)
        """
        return self._element.value_dimension()

    def function_gradient_shape(self):
        """
        :return: shape of the gradient of a shape function (tuple)
        """
        return self._element.value_dimension(), self._element.topological_dimension()

    def evaluate_function(self, function, ref_point):
        return function(self._element.get_mapping().map_point(ref_point))
