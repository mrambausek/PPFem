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

from ppfem.geometry.mesh_entity import MeshEntity, SubEntity


class Line(MeshEntity):

    def __init__(self, vertices, index, mesh):
        MeshEntity.__init__(self, vertices, index, mesh)

    @staticmethod
    def topological_dim():
        return 1

    def set_sub_entities(self):
        """Orientation of vertices is defined to be the axis direction."""
        v1 = self._mesh.select_vertex(self._vertices[0])
        v2 = self._mesh.select_vertex(self._vertices[1])
        if len(v1.coords()) > 1:
            raise Exception("Vertices have more than one coordinate. In such a >1d setting sub entities are not"
                            "defined for lines.")
        if v1[0] < v2[0]:
            self._sub_entities = [SubEntity(v1, SubEntity.neg),
                                  SubEntity(v2, SubEntity.pos)]
        else:
            self._sub_entities = [SubEntity(v1, SubEntity.pos),
                                  SubEntity(v2, SubEntity.neg)]

    def _local_sub_entity_vertices(self):
        # return self._vertices[0], self._vertices[1]
        raise Exception("This method is supposted to be called by 'set_sub_entities' in class 'MeshEntity'."
                        "'set_sub_entities', however, is overriden in class Line and does not call this method.")

