from ppfem.fem.physical_element import PhysicalElement
from ppfem.elements.lagrange_elements import LagrangeLine
from ppfem.geometry.mapping import FEMapping
import scipy as sp


class IsoparametricContinuousLagrange1d(PhysicalElement):

    def __init__(self, dim):
        self._dim = dim
        self._mapping = None
        self._ref_element = None
        self._mesh_entity = None
        self._cache = {}

    def number_of_global_dofs_per_vertex(self):
        return self._dim

    def number_of_global_non_vertex_dofs(self):
        return 0

    def number_of_global_dofs(self):
        return self.number_of_global_dofs_per_vertex() * self._mesh_entity.number_of_vertices() + \
               self.number_of_global_non_vertex_dofs()

    def set_mesh_entity(self, mesh_entity):
        self._mesh_entity = mesh_entity
        vertices = mesh_entity.vertices()
        self._ref_element = LagrangeLine(len(vertices) - 1, self._dim)
        self._mapping = FEMapping(self._ref_element)
        self._mapping.set_mesh_entity(mesh_entity)

    def topological_dimension(self):
        return 1

    def value_dimension(self):
        return self._dim

    def _clear_cache(self):
        self._cache = {}

    def _mapping_is_valid(self, point):
        raise NotImplementedError("This kind of check is not implemented yet.")

    def interpolate_function(self, function):
        return self._ref_element.interpolate_function(function, self._mapping)

    def function_value(self, dof_values, point):
        return self._ref_element.function_value(dof_values, point)

    def function_gradient(self, dof_values, point):
        J_inv = self._mapping.inverse_jacobian(point)
        return self._ref_element.function_gradient(dof_values, point, J_inv)

    def shape_function_values(self, point):
        return self._ref_element.basis_function_values(point)

    def shape_function_gradients(self, point):
        J_inv = self._mapping.inverse_jacobian(point)
        return self._ref_element.basis_function_gradients(point, J_inv)

    def physical_coords(self, ref_point):
        return self._mapping.map_point(ref_point)

    def director(self, point):
        J = self._mapping.jacobian(point)
        return J/sp.linalg.norm(J)

    def jxw(self, qp_data):
        return self._mapping.jacobian_det(qp_data.point) * qp_data.weight

    def global_vertex_indices(self):
        return self._mesh_entity.global_vertex_indices()

    def index(self):
        return self._mesh_entity.index

    def domain_indicator(self):
        return self._mesh_entity.domain_indicator

    def boundary_indicator(self):
        return self._mesh_entity.boundary_indicator
