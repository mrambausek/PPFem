from ppfem.geometry.mesh_entity import MeshEntity


class Face(MeshEntity):

    def __init__(self, vertices, index, mesh):
        MeshEntity.__init__(self, vertices, index, mesh)

    @staticmethod
    def topological_dim():
        return 2
