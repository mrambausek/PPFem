from ppfem.geometry.point import Point
from ppfem.geometry.edge import Edge
from ppfem.geometry.face import Face
from ppfem.geometry.cell import Cell


class Mesh(object):

    def __init__(self, space_dim, topological_dim=None):
        self._cell_dict = {}
        self._face_dict = {}
        self._edge_dict = {}
        self._vertex_dict = {}
        self._space_dim = space_dim
        if topological_dim is None:
            self._topological_dim = space_dim
        else:
            self._topological_dim = topological_dim

    def add_vertex(self, vertex):
        Mesh._add_entity(vertex, vertex.global_index(), self._vertex_dict, "vertex dict")

    def add_face(self, vertex_numbers, number):
        Mesh._add_entity(Face(vertex_numbers, number, self), number, self._face_dict, "face dict")

    def add_edge(self, vertex_numbers, number):
        Mesh._add_entity(Edge(vertex_numbers, number, self), number, self._edge_dict, "edge dict")

    def add_cell(self, vertex_numbers, number):
        Mesh._add_entity(Cell(vertex_numbers, number, self), number, self._cell_dict, "cell dict")

    @staticmethod
    def _add_entity(entity, number, container, container_name):
        if not number in container.keys():
            container[number] = entity
        else:
            raise Exception("There is already an entity with number {0:d} registered in container \"{1:s}\"!"
                            .format(number, container_name))
