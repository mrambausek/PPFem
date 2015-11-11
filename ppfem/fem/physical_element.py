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
    def number_of_global_dofs(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def set_mesh_entity(self, mesh_entity):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def function_value(self, dof_values, point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def function_gradient(self, dof_values, point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def shape_function_values(self, point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def shape_function_gradients(self, dof_values, point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def physical_coords(self, ref_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def jxw(self, quadrature_point):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def global_vertex_indices(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def index(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def domain_indicator(self):
        raise Exception('Abstract method called!')

    @abc.abstractmethod
    def boundary_indicator(self):
        raise Exception('Abstract method called!')
