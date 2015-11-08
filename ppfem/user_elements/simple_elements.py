from ppfem.fem.physical_element import PhysicalElement
from ppfem.elements.lagrange_elements import LagrangeLine
from ppfem.geometry.mapping import FEMapping
from ppfem.geometry.edge import Edge
import scipy as sp


class IsoparametricContinuousLagrange1d(PhysicalElement):

    def __init__(self, material, quadrature, dim):
        self._material = material
        self._quadrature = quadrature
        self._dim = dim
        self._mapping = None
        self._current_entity_index = None
        self._cache = {}

    def number_of_global_dofs_per_vertex(self):
        return self._dim

    def number_of_global_non_vertex_dofs(self):
        return 0

    def _setup_shape(self, mesh_entity):
        vertices = mesh_entity.vertices()
        self._ref_element = LagrangeLine(len(vertices) - 1, self._dim)
        self._mapping = FEMapping(self._ref_element)
        self._mapping.set_mesh_entity(mesh_entity)

    def _clear_cache(self):
        self._cache = {}

    def _mapping_is_valid(self, point):
        raise NotImplementedError("This kind of check is not implemented yet.")

    def _function_value(self, dof_values, point):
        return self._ref_element.function_value(dof_values, point)

    def _function_gradient(self, dof_values, point):
        J_inv = self._mapping.inverse_jacobian(point)
        return self._ref_element.function_gradient(dof_values, point, J_inv)

    def _shape_function_values(self, point):
        return self._ref_element.basis_function_values(point)

    def _shape_function_gradients(self, point):
        J_inv = self._mapping.inverse_jacobian(point)
        return self._ref_element.basis_function_gradients(point, J_inv)

    def _director(self, point):
        J = self._mapping.jacobian(point)
        return J/sp.linalg.norm(J)

    def compute_rhs(self, mesh_entity, dof_values):
        if self._current_entity_index != mesh_entity.index:
            self._clear_cache()
            self._current_entity_index = mesh_entity.index

        if type(mesh_entity) != Edge:
            raise Exception("ContinuousLagrange1d operates on edges!")

        # setup ref. element and mapping
        self._setup_shape(mesh_entity)

        # compute values at quadrature points
        nqp = self._quadrature.number_of_quadrature_points()
        qp_data = self._quadrature.quadrature_data()
        ndofs = self._ref_element.n_of_dofs()
        rhs = sp.zeros(ndofs)

        for qp in qp_data:
            shape_function_values = self._shape_function_values(qp.point)
            shape_function_gradients = self._shape_function_gradients(qp.point)
            function_value = self._function_value(dof_values, qp.point)
            function_gradient = self._function_gradient(dof_values, qp.point)
            physical_coords = self._mapping.map_point(qp.point)
            material_rhs = self._material.compute_rhs(physical_coords, function_value, function_gradient)
            JxW = self._mapping.jacobian_det(qp.point) * qp.weight
            # shape probleme?
            # Suboptimale Struktur des Codes: im Element taucht Wissen über die Linear-Form auf!
            # Das sollte ausgelagert werden!

            # print(shape_function_gradients.shape, shape_function_gradients.T.shape, material_rhs.shape)
            rhs += sp.dot(shape_function_gradients, material_rhs).reshape((ndofs,)) * JxW

        return rhs

    def compute_lhs(self, mesh_entity, dof_values):
        if self._current_entity_index != mesh_entity.index:
            self._clear_cache()
            self._current_entity_index = mesh_entity.index

        if type(mesh_entity) != Edge:
            raise Exception("ContinuousLagrange1d operates on edges!")

        # setup ref. element and mapping
        self._setup_shape(mesh_entity)

        # compute values at quadrature points
        nqp = self._quadrature.number_of_quadrature_points()
        qp_data = self._quadrature.quadrature_data()
        ndofs = self._ref_element.n_of_dofs()
        lhs = sp.zeros((ndofs, ndofs))

        for qp in qp_data:
            shape_function_values = self._shape_function_values(qp.point)
            shape_function_gradients = self._shape_function_gradients(qp.point)
            function_value = self._function_value(dof_values, qp.point)
            function_gradient = self._function_gradient(dof_values, qp.point)
            physical_coords = self._mapping.map_point(qp.point)
            material_lhs = self._material.compute_lhs(physical_coords, function_value, function_gradient)
            JxW = self._mapping.jacobian_det(qp.point) * qp.weight
            # shape probleme?
            # Suboptimale Struktur des Codes: im Element taucht Wissen über die Bilinear-Form auf!
            # Das sollte ausgelagert werden!

            print(shape_function_gradients.shape, material_lhs.shape, shape_function_gradients.T.shape)
            print(shape_function_gradients, material_lhs)
            print(sp.dot(shape_function_gradients, material_lhs), sp.dot(shape_function_gradients, material_lhs).shape)
            print(lhs.shape)
            s = (ndofs, self._dim)
            lhs += sp.dot(sp.dot(shape_function_gradients, material_lhs).reshape(s),
                          shape_function_gradients.reshape(s).T) * JxW

        return lhs
