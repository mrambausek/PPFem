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

from ppfem.geometry.mesh_entity import MeshEntity


class Face(MeshEntity):
    def __init__(self, vertices, index, mesh):
        MeshEntity.__init__(self, vertices, index, mesh)

    @staticmethod
    def topological_dim():
        return 2

    def _local_sub_entity_vertices(self):
        if len(self._vertices) == 3:
            return [(self._vertices[0], self._vertices[1]),
                    (self._vertices[1], self._vertices[2]),
                    (self._vertices[2], self._vertices[0])]
        elif len(self._vertices) == 4:
            return [(self._vertices[0], self._vertices[1]),
                    (self._vertices[1], self._vertices[2]),
                    (self._vertices[2], self._vertices[3]),
                    (self._vertices[3], self._vertices[0])]
        else:
            raise NotImplementedError("Currently only linear triangle and quad faces support setting subentities.")
