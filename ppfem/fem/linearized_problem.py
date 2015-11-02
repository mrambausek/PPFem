import abc


class LinearizedProblem(abs.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def assemble_rhs(self, global_rhs):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_lhs(self, global_lhs):
        raise Exception("Abstract method called!")


class FEAPSytleProblem(LinearizedProblem):

    def __init__(self, mesh, element, boundary_conditions):
        self._mesh = mesh
        self._element = element
        self._boundary_conditions = boundary_conditions
        self._dof_map = {}
        self._global_state_vector = None
        self.number_of_dofs = None

        # TODO:
        # Analyze number of top level mesh entities and global elmt dofs to setup
        # the global dof map (~topology mapping)
        # vertex_index : first index, last index

    def set_global_state(self, new_vector):
        self._global_state_vector = new_vector

    def assemble_rhs(self, global_rhs):
        raise NotImplementedError('Implement me!')

    def assemble_lhs(self, global_lhs):
        raise NotImplementedError('Implement me!')
