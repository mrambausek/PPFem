import abc
import scipy as sp
import scipy.sparse as sparse
from ppfem.fem.dof import DOF


class LinearizedProblem(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def get_sparsity(self, sp_format='coo', mesh_entities=None):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_rhs(self, global_rhs):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_lhs(self, global_lhs):
        raise Exception("Abstract method called!")


class FEAPStyleProblem(LinearizedProblem):
    def __init__(self, mesh, element, boundary_conditions):
        self._mesh = mesh
        self._element = element
        self._boundary_conditions = boundary_conditions
        self._dof_list = {}
        self._element_dof_map = {}
        self._global_state_vector = None
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

        self._global_state_vector = sp.zeros(n_of_dofs)
        self._generate_global_dofs(mesh_entities)
        self.storage_ready = True

    def _get_mesh_entities(self):
        mesh_entities = None
        if self._mesh.topological_dim() == 1:
            mesh_entities = self._mesh.edges()
        elif self._mesh.topological_dim() == 2:
            mesh_entities = self._mesh.faces()
        elif self._mesh.topological_dim() == 3:
            mesh_entities = self._mesh.vertices()
        else:
            raise NotImplementedError("Topological dimension of mesh must be 1, 2 or 3.")
        return mesh_entities

    def _add_dof(self, insert_at, value=0.0, vertex=None):
        self._dof_list[insert_at] = DOF(insert_at, value=value, vertex=vertex)
        return insert_at + 1

    def _add_dofs(self, insert_at, values, vertex=None):
        n = insert_at
        for v in values:
            n = self._add_dof(n, v, vertex)
        return range(insert_at, n), n

    def _generate_global_dofs(self, mesh_entities):
        """Generates and tabulates the global dofs for each element."""
        visited_vertices = []
        insert_at = 0
        vertex_dofs_per_element = self._element.number_of_global_dofs_per_vertex()
        non_vertex_dofs_per_element = self._element.number_of_global_non_vertex_dofs()

        vertex_zeros = sp.zeros(vertex_dofs_per_element)
        non_vertex_zeros = sp.zeros(non_vertex_dofs_per_element)

        for e in mesh_entities:
            dofs = []
            for v in e.global_vertex_indices():
                if v in visited_vertices:
                    continue
                visited_vertices.append(v)
                new_dofs, insert_at = self._add_dofs(insert_at, vertex_zeros, v)
                self._mesh.vertex(v).set_associated_dof_numbers(new_dofs)
                dofs += new_dofs

            new_dofs, insert_at = self._add_dofs(insert_at, non_vertex_zeros)
            dofs += new_dofs
            self._element_dof_map[e.index] = dofs

    def _get_element_dof_index_array(self, element):
        return sp.array(self._element_dof_map[element.index], dtype=sp.int64)

    def _get_element_dof_values_array(self, element):
        return sp.array([self._global_state_vector[dof] for dof in self._element_dof_map[element.index]])

    def _get_element_coo_assembly_indices(self, element, global_data=None):
        e_dofs = self._get_element_dof_index_array(element)
        m = sp.vstack([e_dofs] * e_dofs.size)
        if global_data is None:
            rows_add, cols_add = m.T.flatten().tolist(), m.flatten().tolist()
        else:
            global_coo_row_indices, global_coo_col_indices = global_data
            rows_add, cols_add = [], []
            row_candidates, col_candidates = m.T.flatten().tolist(), m.flatten().tolist()

            for i in range(len(global_coo_row_indices)):
                p_global = [global_coo_row_indices[i], global_coo_col_indices[i]]
                for p_local in zip(row_candidates, col_candidates):
                    if p_global[0] != p_local[0] or p_global[1] != p_local[1]:
                        rows_add.append(p_local[0])
                        cols_add.append(p_local[1])
        return rows_add, cols_add

    def _get_sparsity_coo(self, mesh_entities, condense=False):
        if mesh_entities is None:
            mesh_entities = self._get_mesh_entities()
        global_coo_row_indices = []
        global_coo_col_indices = []

        for e in mesh_entities:
            if condense:
                new_pairs = self._get_element_coo_assembly_indices(e, global_data=(global_coo_row_indices,
                                                                                   global_coo_col_indices))
            else:
                new_pairs = self._get_element_coo_assembly_indices(e)

            global_coo_row_indices += new_pairs[0]
            global_coo_col_indices += new_pairs[1]
        return sp.array(global_coo_row_indices, dtype=sp.int64), sp.array(global_coo_col_indices, dtype=sp.int64)

    def set_global_state(self, new_vector):
        self._global_state_vector = new_vector

    def _get_element_rhs(self, element):
        dof_values = self._get_element_dof_values_array(element)
        return self._element.compute_rhs(element, dof_values)

    def _get_element_matrix(self, element):
        dof_values = self._get_element_dof_values_array(element)
        return self._element.compute_lhs(element, dof_values)

    def _assemble_local_rhs(self, element, global_rhs):
        k = self._get_element_matrix(element)
        dofs = self._get_element_dof_index_array(element)
        i = 0
        for I in dofs:
            global_rhs[I] = k[i]
            i += 1

    def _assemble_local_matrix(self, element, global_lhs):
        k = self._get_element_matrix(element)
        dofs = self._get_element_dof_index_array(element)
        i, j = 0, 0
        for I in dofs:
            for J in dofs:
                global_lhs[I,J] = k[i,j]
                j += 1
            i += 1

    def get_sparsity(self, sp_format='coo', mesh_entities=None):
        if not self.storage_ready:
            raise Exception('Storage/DoFs not ready yet!')

        if mesh_entities is None:
            mesh_entities = self._get_mesh_entities()

        if sp_format == 'coo':
            return self._get_sparsity_coo(mesh_entities)
        elif sp_format == 'coo_condensed':
            return self._get_sparsity_coo(mesh_entities, condense=True)

    def assemble_rhs(self, global_rhs):
        for e in self._get_mesh_entities():
            self._assemble_local_rhs(e, global_rhs)

    def assemble_lhs(self, global_lhs):
        for e in self._get_mesh_entities():
            self._assemble_local_matrix(e, global_lhs)
