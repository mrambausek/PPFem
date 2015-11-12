import scipy as sp


class FunctionSpace(object):

    def __init__(self, mesh, element_map):
        self._element_map = element_map
        self._mesh = mesh
        self._element_dof_map = {}
        self._vertex_dof_map = {}
        self.number_of_dofs = None

        self.storage_ready = False
        self._setup_storage()

    def _setup_storage(self):
        for v in self._mesh.vertices():
            self._vertex_dof_map[v.global_index()] = []

        self.number_of_dofs = self._generate_assembly_data()
        self.storage_ready = True

    def get_mesh_entities(self):
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

    def _generate_element_vertex_dofs(self, element, first_dof_index, visited_vertices):
        vertex_dofs_per_element = element.number_of_global_dofs_per_vertex()

        dofs = []
        for v in element.global_vertex_indices():
            if v not in visited_vertices:
                next_first_dof_index = first_dof_index + vertex_dofs_per_element
                new_dofs = range(first_dof_index, next_first_dof_index)
                first_dof_index = next_first_dof_index
                self._vertex_dof_map[v] += new_dofs
            dofs += self._vertex_dof_map[v]

        return dofs, first_dof_index

    def _generate_element_non_vertex_dofs(self, element, first_dof_index):
        non_vertex_dofs_per_element = element.number_of_global_non_vertex_dofs()
        next_first_dof_index = first_dof_index+non_vertex_dofs_per_element
        dofs = range(first_dof_index, first_dof_index+non_vertex_dofs_per_element)
        return dofs, next_first_dof_index

    def _generate_assembly_data(self):
        """Generates and tabulates the global dofs for each element."""
        visited_vertices = []
        first_dof_index = 0

        for e in self.get_mesh_entities():
            element = self.get_element(e)

            dofs = []
            new_dofs, first_dof_index = self._generate_element_vertex_dofs(
                element, first_dof_index, visited_vertices
            )
            dofs += new_dofs
            visited_vertices += e.global_vertex_indices()

            new_dofs, first_dof_index = self._generate_element_non_vertex_dofs(element, first_dof_index)
            dofs += new_dofs

            self._element_dof_map[element.index()] = dofs
        return first_dof_index

    def get_element_dof_index_array(self, element_index):
        return sp.array(self._element_dof_map[element_index], dtype=sp.int64)

    def get_element(self, mesh_entity):
        elmt = self._element_map[mesh_entity.domain_indicator]
        elmt.set_mesh_entity(mesh_entity)
        return elmt

    def function_dim(self):
        pass

    def function_space_dim(self):
        pass


class TrialFunctionSpace(FunctionSpace):
    def __init__(self, mesh, element_map):
        FunctionSpace.__init__(mesh, element_map)
        self._global_dof_vector = sp.zeros(self.number_of_dofs)

    def set_global_dofs(self, new_vector):
        self._global_dof_vector = new_vector

    def get_element_dof_values_array(self, element_index):
        self._global_dof_vector[self._element_dof_map[element_index]]
