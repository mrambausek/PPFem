from ppfem.geometry.mesh_entity import MeshEntity


class Line(MeshEntity):

    def __init__(self, vertices, index, mesh):
        MeshEntity.__init__(self, vertices, index, mesh)

    @staticmethod
    def topological_dim():
        return 1
