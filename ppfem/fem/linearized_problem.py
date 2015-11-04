import abc
import scipy as sp
import scipy.sparse as sparse
from ppfem.fem.dof import DOF


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
        self._dof_list = {}
        self._element_dof_map = {}
        self._global_state_vector = None
        self._global_matrix = None
        self._global_rhs = None
        self.number_of_dofs = None

        # for coo-matrix format: fast but high storage demands!
        self._global_coo_row_indices = None
        self._global_coo_col_indices = None
        self._global_coo_matrix_data = None

        self._matrix_assembly_format = 'coo'
        self.storage_ready = False
        self.setup_storage()

    def setup_storage(self):
        mesh_entities = self._get_mesh_entities()
        n_of_dofs = len(self._mesh.vertices()) * self._element.number_of_global_dofs_per_vertex() \
                  + len(mesh_entities) * self._element.number_of_global_non_vertex_dofs()
        self._dof_list = [None] * n_of_dofs

        self._global_rhs = sp.zeros(n_of_dofs)
        self._global_state_vector = sp.zeros(n_of_dofs)

        self._generate_global_dofs(mesh_entities)

        if self._matrix_assembly_format == 'coo':
            self._generate_coo_lhs_sparsity(mesh_entities)
        else:
            raise NotImplementedError("(Pre-)Assembly format {0:s} is not supported."
                                      .format(self._matrix_assembly_format))

        self.storage_ready = True

    def _get_mesh_entities(self):
        mesh_entities = None
        if self._mesh.topological_dim() == 1:
            mesh_entities = self._mesh.edges()
        elif self._mesh.topological_dim() == 2:
            mesh_entities = self._mesh.faces()
        elif self._mesh.topological_dim() == 3:
            mesh_entities = self._mesh.vertives()
        else:
            raise NotImplementedError("Topological dimension of mesh must be 1, 2 or 3.")
        return mesh_entities

    def _add_dof(self, number, value=0.0, vertex=None):
        self._dof_list[number] = DOF(number, value=value, vertex=vertex)

    def _add_dofs(self, first_number, values, vertex=None):
        n = first_number
        for v in values:
            self._add_dof(n, v, vertex)
            n += 1
        return range(first_number, n+1)

    def _generate_global_dofs(self, mesh_entities):
        """Generates and tabulates the global dofs for each element."""
        visited_vertices = []
        current_number = 0
        vertex_dofs_per_element = self._element.number_of_global_dofs_per_vertex()
        non_vertex_dofs_per_element = self._element.number_of_global_non_vertex_dofs()

        vertex_zeros = sp.zeros(vertex_dofs_per_element)
        non_vertex_zeros = sp.zeros(non_vertex_dofs_per_element)

        for e in mesh_entities:
            dofs = []
            for v in e.global_vertex_indices:
                if v in visited_vertices:
                    continue
                visited_vertices.append(v)
                dofs += self._add_dofs(current_number, vertex_zeros, v)
                current_number = dofs[-1]

            dofs += self._add_dofs(current_number, non_vertex_zeros)
            current_number = dofs[-1]
            self._element_dof_map[e.index] = dofs

    def _generate_lhs_coo_sparsity(self, mesh_entitites):
        self._global_coo_row_indices = []
        self._global_coo_col_indices = []

        for e in mesh_entitites:
            e_dofs = sp.array( self._element_dof_map[e.index], dtype=sp.int64 )
            m = sp.vstack( [e_dofs]*e_dofs.size )
            self._global_coo_row_indices += m.T.flatten().tolist()
            self._global_coo_col_indices += m.flatten().tolist()

        self._global_coo_row_indices = sp.array(self._global_coo_row_indices, dtype=sp.int64)
        self._global_coo_col_indices = sp.array(self._global_coo_col_indices, dtype=sp.int64)

    def set_global_state(self, new_vector):
        self._global_state_vector = new_vector

    def assemble_rhs(self, global_rhs):
        raise NotImplementedError('Implement me!')

    def assemble_lhs(self, global_lhs):
        raise NotImplementedError('Implement me!')
