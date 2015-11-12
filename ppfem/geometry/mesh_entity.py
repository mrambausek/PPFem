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
import scipy as sp


class MeshEntity(abc.ABC):

    def __init__(self, vertices, index, mesh):
        self._vertices = tuple(vertices)
        self._number_of_vertices = len(self._vertices)
        self.index = index
        self.domain_indicator = 0
        self.boundary_indicator = None
        self._mesh = mesh

    def global_vertex_indices(self):
        return self._vertices

    def number_of_vertices(self):
        return self._number_of_vertices

    def vertices(self):
        return self._mesh.select_vertices(self._vertices)

    def vertex_coords(self, local_vertex_number=None):
        if local_vertex_number is None:
            return sp.array([v.coords() for v in self.vertices()])
        else:
            return self._mesh.select_vertex(local_vertex_number).coords()

    def vertex(self, local_vertex_number):
        return self._mesh.vertex(self._vertices[local_vertex_number])

    @staticmethod
    @abc.abstractmethod
    def topological_dim():
        raise Exception('Abstract method called!')
