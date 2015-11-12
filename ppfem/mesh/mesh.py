from ppfem.geometry.point import Point
from ppfem.geometry.line import Line
from ppfem.geometry.face import Face
from ppfem.geometry.cell import Cell


class Mesh(object):

    def __init__(self, space_dim, topological_dim=None):
        self._cell_dict = {}
        self._face_dict = {}
        self._line_elem_dict = {}
        self._vertex_dict = {}
        self._space_dim = space_dim
        if topological_dim is None:
            self._topological_dim = space_dim
        else:
            self._topological_dim = topological_dim

    def get_mesh_entities(self):
        if self.topological_dim() == 1:
            return self.lines()
        elif self.topological_dim() == 2:
            return self.faces()
        elif self.topological_dim() == 3:
            return self.vertices()
        else:
            raise NotImplementedError("Topological dimension of mesh must be 1, 2 or 3.")

    def add_vertex(self, vertex):
        Mesh._add_entity(vertex, vertex.global_index(), self._vertex_dict, "vertex dict")

    def add_face(self, vertex_numbers, number):
        Mesh._add_entity(Face(vertex_numbers, number, self), number, self._face_dict, "face dict")

    def add_line(self, vertex_numbers, number):
        Mesh._add_entity(Line(vertex_numbers, number, self), number, self._line_elem_dict, "edge dict")

    def add_cell(self, vertex_numbers, number):
        Mesh._add_entity(Cell(vertex_numbers, number, self), number, self._cell_dict, "cell dict")

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
        return self._line_elem_dict.values()

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
