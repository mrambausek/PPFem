import scipy as sp


class FunctionSpace(object):

    def __init__(self, element, mesh=None, subdomain=None):
        self._element = element
        self._mesh = mesh
        self._subdomain = subdomain
        self._element_dof_map = {}
        self._vertex_dof_map = {}
        self.number_of_dofs = None

        self.storage_ready = False
        if mesh is not None:
            self._setup_storage()

    def _setup_storage(self):
        for v in self._mesh.vertices():
            self._vertex_dof_map[v.global_index()] = []

        self.number_of_dofs = self._generate_assembly_data()
        self.storage_ready = True

    def set_mesh(self, mesh, sub_domain=None):
        self._mesh = mesh
        self._subdomain = sub_domain
        self._setup_storage()

    def mesh_entity_iterator(self):
        if self._subdomain is None:
            return self._mesh.get_mesh_entities()
        return filter(lambda e: e.domain_indicator == self._subdomain, self._mesh.get_mesh_entities())

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

        for e in self._mesh.get_mesh_entities():
            if self._subdomain is not None and e.domain_indicator != self._subdomain:
                continue

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
        if self._subdomain is not None and mesh_entity.domain_indicator != self._subdomain:
            raise Exception("Mesh entity in wrong subdomain!")
        self._element.set_mesh_entity(mesh_entity)
        return self._element

    def localize(self, mesh_entity):
        return LocalFunctionSpace(self.get_element(mesh_entity))

    def function_dim(self):
        self._element.value_dim()

    def function_space_dim(self):
        return self.number_of_dofs

    def interpolate_function(self, function):
        global_dofs = sp.zeros(self.number_of_dofs)
        for e in self.mesh_entity_iterator():
            global_dofs[self._element_dof_map[e.index]] = self.localize(e).interpolate_function(function)
        return global_dofs


class LocalFunctionSpace(object):
    def __init__(self, element):
        self._element = element

    def interpolate_function(self, function):
        return self._element.interpolate_funtion(function)

    def shape_function_values(self, point):
        self._element.shape_function_values(point)

    def shape_function_gradients(self, point):
        self._element.shape_function_gradients(point)

    def physical_coords(self, ref_point):
        self._element.phyical_coords(ref_point)

    def jxw(self, qp_data):
        self._element.jxw(qp_data)
