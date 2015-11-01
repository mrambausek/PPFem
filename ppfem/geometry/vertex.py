from ppfem.geometry.point import Point


class Vertex(Point):

    def __init__(self, *coords, global_number):
        Point.__init__(self, coords, global_number)

    def global_index(self):
        return self.index

    def is_part_of(self, mesh_entitiy):
        return self.index in mesh_entitiy.global_vertex_indices()
