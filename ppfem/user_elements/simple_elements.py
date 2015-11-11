from ppfem.fem.physical_element import PhysicalElement
from ppfem.elements.lagrange_elements import LagrangeLine
from ppfem.geometry.mapping import FEMapping
import scipy as sp


class LazyEval(object):
    def __init__(self, func, args):
        self._func = func
        self._args = args

    def __call__(self, *args, **kwargs):
        return self._func(*self._args)


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

    def _clear_cache(self):
        self._cache = {}

    def _mapping_is_valid(self, point):
        raise NotImplementedError("This kind of check is not implemented yet.")

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

    def jxw(self, qp):
        return self._mapping.jacobian_det(qp.point) * qp.weight

    def global_vertex_indices(self):
        return self._mesh_entity.global_vertex_indices()

    def index(self):
        return self._mesh_entity.index

    def domain_indicator(self):
        return self._mesh_entity.domain_indicator

    def boundary_indicator(self):
        return self._mesh_entity.boundary_indicator

    # def _point_data(self, dof_values, q_point):
    #     return EvalData({
    #         'physical_coords': LazyEval(self._mapping.map_point, (q_point,)),
    #         'function_value': LazyEval(self._function_value, (dof_values, q_point)),
    #         'function_gradient': LazyEval(self._function_gradient, (dof_values, q_point))
    #     })

    # def compute_rhs(self, mesh_entity, dof_values):
    #     if self._current_entity_index != mesh_entity.index:
    #         self._clear_cache()
    #         self._current_entity_index = mesh_entity.index
    #
    #     if type(mesh_entity) != Edge:
    #         raise Exception("ContinuousLagrange1d operates on edges!")
    #
    #     # setup ref. element and mapping
    #     self._setup_shape(mesh_entity)
    #
    #     # compute values at quadrature points
    #     ndofs = self._ref_element.n_of_dofs()
    #     rhs = sp.zeros(ndofs)
    #
    #     rhs_contraction = self._material.contraction_data().rhs
    #
    #     # quadrature loop
    #     for qp in self._quadrature.quadrature_data():
    #         material_rhs = self._material.compute_rhs(self._point_data(dof_values, qp.point))
    #
    #         if rhs_contraction == ContractionData.func:
    #             shape_function_values = self._shape_function_values(qp.point)
    #             rhs += \
    #                 sp.dot(shape_function_values, material_rhs)\
    #                 .reshape((ndofs,)) \
    #                 * self._jxw(qp)
    #
    #         elif rhs_contraction == ContractionData.grad:
    #             shape_function_gradients = self._shape_function_gradients(qp.point)
    #             rhs += \
    #                 sp.dot(shape_function_gradients, material_rhs)\
    #                 .reshape((ndofs,)) \
    #                 * self._jxw(qp)
    #         else:
    #             raise NotImplementedError("rhs_contraction \"{:s} is not implemented.\""
    #                                       .format(rhs_contraction))
    #
    #     return rhs
    #
    # def compute_lhs(self, mesh_entity, dof_values):
    #     if self._current_entity_index != mesh_entity.index:
    #         self._clear_cache()
    #         self._current_entity_index = mesh_entity.index
    #
    #     if type(mesh_entity) != Edge:
    #         raise Exception("ContinuousLagrange1d operates on edges!")
    #
    #     # setup ref. element and mapping
    #     self._setup_shape(mesh_entity)
    #
    #     # compute values at quadrature points
    #     ndofs = self._ref_element.n_of_dofs()
    #     shape = (ndofs, self._dim)  # in 1d, shape of value and gradient is the same!
    #     lhs = sp.zeros((ndofs, ndofs))
    #
    #     lhs_contraction = self._material.contraction_data().lhs
    #
    #     # quadrature loop
    #     for qp in self._quadrature.quadrature_data():
    #         material_lhs = self._material.compute_lhs(self._point_data(dof_values, qp.point))
    #
    #         if lhs_contraction == ContractionData.func:
    #             shape_function_values = self._shape_function_values(qp.point)
    #             lhs += \
    #                 sp.dot(sp.dot(shape_function_values, material_lhs)
    #                        .reshape(shape),
    #                        shape_function_values.reshape(shape).T) \
    #                 * self._jxw(qp)
    #
    #         elif lhs_contraction == ContractionData.grad:
    #             shape_function_gradients = self._shape_function_gradients(qp.point)
    #             lhs += \
    #                 sp.dot(sp.dot(shape_function_gradients, material_lhs)
    #                        .reshape(shape),
    #                        shape_function_gradients.reshape(shape).T) \
    #                 * self._jxw(qp)
    #         else:
    #             raise NotImplementedError("lhs_contraction \"{:s} is not implemented.\""
    #                                       .format(lhs_contraction))
    #
    #     return lhs
