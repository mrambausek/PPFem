from ppfem.fem.assembler import Assembler, find_index_pair
import scipy as sp
from ppfem.fem.dof import DOF


class FiniteElementModel(Assembler):
    smf_coo = 0
    smf_coo_condensed = 1

    def __init__(self, mesh, dgl_map, boundary_conditions):
        self._mesh = mesh
        self._dgl_map = dgl_map
        self._boundary_conditions = boundary_conditions
        self._dof_list = []
        self._element_dof_map = {}
        self._vertex_dof_map = {}
        self._global_state_vector = None
        self.number_of_dofs = None

        self.storage_ready = False
        self.setup_storage()

    def setup_storage(self):
        for v in self._mesh.vertices():
            self._vertex_dof_map[v.global_index()] = []

        self.number_of_dofs = self._generate_global_dofs()
        self._global_state_vector = sp.zeros(self.number_of_dofs)
        self.storage_ready = True

    def _get_mesh_entities(self):
        mesh_entities = None
        if self._mesh.topological_dim() == 1:
            mesh_entities = self._mesh.lines()
        elif self._mesh.topological_dim() == 2:
            mesh_entities = self._mesh.faces()
        elif self._mesh.topological_dim() == 3:
            mesh_entities = self._mesh.vertices()
        else:
            raise NotImplementedError("Topological dimension of mesh must be 1, 2 or 3.")
        return mesh_entities

    def _add_dof(self, insert_at, value=0.0, vertex=None, associate_to_vertex=True):
        self._dof_list[insert_at] = DOF(insert_at, value=value, vertex=vertex)
        if vertex is not None:
            self._vertex_dof_map[vertex].append(insert_at)
        return insert_at + 1

    def _add_dofs(self, insert_at, values, vertex=None):
        n = insert_at
        for v in values:
            n = self._add_dof(n, v)
        new_dofs = range(insert_at, n)
        if vertex is not None:
            self._vertex_dof_map[vertex] += new_dofs
        return new_dofs, n

    def _generate_element_vertex_dofs(self, element, insert_at, visited_vertices):
        vertex_dofs_per_element = element.number_of_global_dofs_per_vertex()

        dofs = []
        for v in element.global_vertex_indices():
            if v in visited_vertices:
                dofs += self._vertex_dof_map[v]
                continue
            new_dofs, insert_at = self._add_dofs(insert_at,
                                                 sp.zeros(vertex_dofs_per_element),
                                                 vertex=v)
            dofs += new_dofs
        return dofs, insert_at

    def _generate_element_non_vertex_dofs(self, element, insert_at):
        non_vertex_dofs_per_element = element.number_of_global_non_vertex_dofs()
        return self._add_dofs(insert_at, sp.zeros(non_vertex_dofs_per_element))

    def _generate_global_dofs(self):
        """Generates and tabulates the global dofs for each element."""
        visited_vertices = []
        insert_at = 0

        for e in self._get_mesh_entities():
            element = self._get_element(e)

            self._dof_list += [None] * element.number_of_global_dofs()
            dofs = []
            print(insert_at, len(self._dof_list))
            new_dofs, insert_at = self._generate_element_vertex_dofs(
                element, insert_at, visited_vertices
            )
            dofs += new_dofs
            visited_vertices += e.global_vertex_indices()

            new_dofs, insert_at = self._generate_element_non_vertex_dofs(element, insert_at)
            dofs += new_dofs

            self._element_dof_map[element.index()] = dofs
        return len(self._dof_list)

    def _get_element_dof_index_array(self, element_index):
        return sp.array(self._element_dof_map[element_index], dtype=sp.int64)

    def _get_element_dof_values_array(self, element_index):
        return sp.array([self._global_state_vector[dof] for dof in self._element_dof_map[element_index]])

    def _get_element_coo_assembly_indices(self, element_index, global_data=None):
        e_dofs = self._get_element_dof_index_array(element_index)
        m = sp.vstack([e_dofs] * e_dofs.size)
        if global_data is None:
            rows_add, cols_add = m.T.flatten().tolist(), m.flatten().tolist()
        else:
            rows_add, cols_add = [], []
            row_candidates, col_candidates = m.T.flatten().tolist(), m.flatten().tolist()
            for p_local in zip(row_candidates, col_candidates):
                if not find_index_pair(p_local, global_data):
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
                new_pairs = self._get_element_coo_assembly_indices(e.index,
                                                                   global_data=(global_coo_row_indices,
                                                                                global_coo_col_indices))
            else:
                new_pairs = self._get_element_coo_assembly_indices(e.index)

            global_coo_row_indices += new_pairs[0]
            global_coo_col_indices += new_pairs[1]
        return sp.array(global_coo_row_indices, dtype=sp.int64), sp.array(global_coo_col_indices, dtype=sp.int64)

    def _get_element(self, mesh_entity):
        return self._dgl_map[mesh_entity.domain_indicator].element(mesh_entity)

    def set_global_state(self, new_vector):
        self._global_state_vector = new_vector

    def _extract_element_dof_values(self, element_index):
        return self._global_state_vector[self._element_dof_map[element_index]]

    def _assemble_local_rhs(self, mesh_entity, global_rhs):
        di = mesh_entity.domain_indicator
        dgl = self._dgl_map[di]
        dof_values = self._get_element_dof_values_array(mesh_entity.index)

        r = dgl.linear_form(mesh_entity, dof_values, None)
        dofs = self._get_element_dof_index_array(mesh_entity.index)
        i = 0
        for I in dofs:
            global_rhs[I] = r[i]
            i += 1

    def _assemble_local_matrix(self, mesh_entity, global_lhs):
        di = mesh_entity.domain_indicator
        dgl = self._dgl_map[di]
        dof_values = self._get_element_dof_values_array(mesh_entity.index)

        k = dgl.bilinear_form(mesh_entity, dof_values, None)
        dofs = self._get_element_dof_index_array(mesh_entity.index)

        i = 0
        for I in dofs:
            j = 0
            for J in dofs:
                global_lhs[I,J] = k[i,j]
                j += 1
            i += 1

    def get_sparsity(self, sp_format=0, mesh_entities=None):
        if not self.storage_ready:
            raise Exception('Storage/DoFs not ready yet!')

        if mesh_entities is None:
            mesh_entities = self._get_mesh_entities()

        if sp_format == FiniteElementModel.smf_coo:
            return self._get_sparsity_coo(mesh_entities)
        elif sp_format == FiniteElementModel.smf_coo_condensed:
            return self._get_sparsity_coo(mesh_entities, condense=True)

    def assemble_rhs(self, global_rhs):
        for e in self._get_mesh_entities():
            self._assemble_local_rhs(e, global_rhs)

    def assemble_lhs(self, global_lhs):
        for e in self._get_mesh_entities():
            self._assemble_local_matrix(e, global_lhs)
