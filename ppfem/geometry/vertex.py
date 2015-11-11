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

    # def associate_dof_number(self, dof_number):
    #     self._associated_dof_numbers.append(dof_number)
    #
    # def set_associated_dof_numbers(self, global_dof_numbers):
    #     self._associated_dof_numbers = global_dof_numbers
    #
    # def associated_dof_numbers(self):
    #     return self._associated_dof_numbers
