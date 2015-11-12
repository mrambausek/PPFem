import abc


class PDE(abc.ABC):
    def __init__(self, element, quadrature):
        self._element = element
        self._quadrature = quadrature

    def element(self, mesh_entity=None):
        if mesh_entity is not None:
            self._element.set_mesh_entity(mesh_entity)
        return self._element

    def quadrature(self):
        return self._quadrature

    def set_element(self, element):
        self._element = element

    def set_quadrature(self, quadrature):
        self._quadrature = quadrature

    @abc.abstractmethod
    def linear_form(self, mesh_entity, dof_values, params):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def bilinear_form(self, mesh_entity, dof_values, params):
        raise Exception("Abstract method called!")
