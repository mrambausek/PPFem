import abc


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
        raise NotImplementedError("Needs to be implemented!")

    def vertex(self, local_vertex_number):
        return self._mesh.vertex(self._vertices[local_vertex_number])

    @staticmethod
    @abc.abstractmethod
    def topological_dim():
        raise Exception('Abstract method called!')
