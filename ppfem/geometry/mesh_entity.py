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
import scipy as sp


class MeshEntity(abc.ABC):
    def __init__(self, vertices, index, mesh):
        self._vertices = tuple(vertices)
        self._number_of_vertices = len(self._vertices)
        self._mesh = mesh
        self._sub_entities = None
        self.index = index
        self.domain_indicator = 0
        self.boundary_indicator = None
        # self.mapping_index = None
        # self.quadrature_index = None

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

    @abc.abstractmethod
    def _local_sub_entity_vertices(self):
        raise Exception('Abstract method called!')

    @staticmethod
    @abc.abstractmethod
    def topological_dim():
        raise Exception('Abstract method called!')

    def orientation(self, vertex_indices):
        pos1 = self._vertices.index(vertex_indices[0])
        pos2 = self._vertices.index(vertex_indices[1])
        if pos1 > pos2:
            return SubEntity.neg
        else:
            return SubEntity.pos

    def set_sub_entities(self):
        self._sub_entities = [self._get_sub_data(local_vertices)
                              for local_vertices in self._local_sub_entity_vertices()]

    def _sub_entity_indices(self):
        return [s.index for s in self._sub_entities]

    def sub_entities(self, filter_func=None):
        if self._sub_entities is None:
            return None
        elif filter_func is None:
            return self._sub_entities
        else:
            return filter(filter_func, self._sub_entities)
        # subs = sorted(self._mesh.get_mesh_entities(topological_dim=self.topological_dim() - 1,
        #                                            indices=self._sub_entity_indices(),
        #                                            filter=filter_func),
        #               key=self._sub_entity_indices().index)

    def get_sub_entity(self, sub_entity_index):
        return self._sub_entities[sub_entity_index]

    def _get_sub_data(self, local_sub_entity_vertices):
        found = list(self._mesh.find_entities_with_vertices(local_sub_entity_vertices, self.topological_dim() - 1))
        if len(found) == 0:
            found.append(self._mesh.add_line(local_sub_entity_vertices))
        elif len(found) > 1:
            raise Exception("More than one entity with the given vertex indices were found!")

        sub_data = SubEntity(found[0], found[0].ortientation(local_sub_entity_vertices))
        return sub_data


class SubEntity(object):
    pos = 1
    neg = -1

    def __init__(self, entity, orientation):
        # entity not necessarily a 'real' mesh entity (Line, Face, Cell); could be a vertex too
        self.entity = entity
        self.orientation = orientation
