import abc
import scipy as sp


class MeshEntity(abc.ABC):

    def __init__(self, vertices, index, mesh):
        self._vertices = tuple(vertices)
        self.index = index
        self._mesh = mesh

    def global_vertex_indices(self):
        return self._vertices

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
