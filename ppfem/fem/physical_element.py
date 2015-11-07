import abc


class PhysicalElement(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def number_of_global_dofs_per_vertex(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def number_of_global_non_vertex_dofs(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def compute_rhs(self, mesh_entity, dof_values):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def compute_lhs(self, mesh_entity, dof_values):
        raise Exception('Abstract method called!')
