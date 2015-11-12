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

from ppfem.geometry.point import Point


class Vertex(Point):

    def __init__(self, coords, global_number, associated_dof_numbers=None):
        Point.__init__(self, *coords, index=global_number)
        if associated_dof_numbers is None:
            self._associated_dof_numbers = []
        else:
            self._associated_dof_numbers = associated_dof_numbers

    def global_index(self):
        return self.index

    def is_part_of(self, mesh_entitiy):
        return self.index in mesh_entitiy.global_vertex_indices()
