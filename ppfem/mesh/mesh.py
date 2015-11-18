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

from ppfem.geometry.line import Line
from ppfem.geometry.face import Face
from ppfem.geometry.cell import Cell


class FilterIndices(object):
    def __init__(self, indices):
        self.indices = indices

    def __call__(self, e):
        return e.index in self.indices


class FilterEntitiesOnVertexIndices(object):
    def __init__(self, vertex_indices):
        self.vertex_indices = frozenset(vertex_indices)

    def __call__(self, e):
        self.vertex_indices == frozenset(e.global_vertex_indices())


class Mesh(object):

    def __init__(self, space_dim, topological_dim=None):
        self._cell_dict = {}
        self._face_dict = {}
        self._line_dict = {}
        self._vertex_dict = {}
        self._space_dim = space_dim
        if topological_dim is None:
            self._topological_dim = space_dim
        else:
            self._topological_dim = topological_dim

    def get_mesh_entities(self, topological_dim=None, indices=None, filter_func=None):
        _iters = []
        if topological_dim is None:
            topological_dim = self.topological_dim()
        if topological_dim == 0:
            _iter = self.vertices()
        elif topological_dim == 1:
            _iter = self.lines()
        elif topological_dim == 2:
            _iter = self.faces()
        elif topological_dim == 3:
            _iter = self.cells()
        else:
            raise NotImplementedError("Topological dimension of mesh entities must be 0, 1, 2 or 3.")

        _iters.append(_iter)

        if filter_func is not None:
            _iters.append(filter(filter_func, _iters[-1]))
        if indices is not None:
            _iters.append(filter(FilterIndices(indices), _iters[-1]))
        return _iters[-1]

    def add_vertex(self, vertex):
        if vertex.index is None:
            if len(self._vertex_dict) == 0:
                vertex.index = 0
            else:
                vertex.index = max(self._vertex_dict.keys()) + 1
        Mesh._add_entity(vertex, vertex.global_index(), self._vertex_dict, "vertex dict")
        return vertex

    def add_line(self, vertex_numbers, number=None):
        if number is None:
            if len(self._line_dict) == 0:
                number = 0
            else:
                number = max(self._line_dict.keys()) + 1
        Mesh._add_entity(Line(vertex_numbers, number, self), number, self._line_dict, "edge dict")
        return self._line_dict[number]

    def add_face(self, vertex_numbers, number=None):
        if number is None:
            if len(self._line_dict) == 0:
                number = 0
            else:
                number = max(self._face_dict.keys()) + 1
        Mesh._add_entity(Face(vertex_numbers, number, self), number, self._face_dict, "face dict")
        return self._face_dict[number]

    def add_cell(self, vertex_numbers, number=None):
        if number is None:
            if len(self._line_dict) == 0:
                number = 0
            else:
                number = max(self._line_dict.keys()) + 1
        Mesh._add_entity(Cell(vertex_numbers, number, self), number, self._cell_dict, "cell dict")
        return self._cell_dict[number]

    def find_entities_with_vertices(self, vertex_indices, topological_dim):
        if topological_dim == 1:
            return filter(FilterEntitiesOnVertexIndices(vertex_indices), self.lines())[0]
        elif topological_dim == 2:
            return filter(FilterEntitiesOnVertexIndices(vertex_indices), self.faces())
        elif topological_dim == 3:
            return filter(FilterEntitiesOnVertexIndices(vertex_indices), self.cell())
        else:
            raise NotImplementedError("Topological dimension of mesh entities must be 1, 2 or 3.")

    def select_vertex(self, global_vertex_number):
        return self._vertex_dict[global_vertex_number]

    def select_vertices(self, global_vertex_numbers):
        return [self._vertex_dict[v] for v in global_vertex_numbers]

    def vertex(self, global_vertex_number):
        return self._vertex_dict[global_vertex_number]

    def vertices(self):
        return self._vertex_dict.values()

    def cells(self):
        return self._cell_dict.values()

    def faces(self):
        return self._face_dict.values()

    def lines(self):
        return self._line_dict.values()

    def topological_dim(self):
        return self._topological_dim

    def space_dim(self):
        return self._space_dim

    @staticmethod
    def _add_entity(entity, number, container, container_name):
        if not number in container.keys():
            container[number] = entity
        else:
            raise Exception("There is already an entity with number {0:d} registered in container \"{1:s}\"!"
                            .format(number, container_name))
