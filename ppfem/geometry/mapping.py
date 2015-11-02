import abc
import numpy as np
import numpy.linalg as npl


class Mapping(abc.ABC):

    def __init__(self, reference_element):
        self._reference_element = reference_element
        self._mesh_entity = None

    def set_mesh_entity(self, mesh_entity):
        self._mesh_entity

    def clear_mesh_entity(self):
        self._mesh_entity = None

    @abc.abstractmethod
    def map_point(self, reference_point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def map_points(self, reference_points):
        return [self.map_point(p) for p in reference_points]

    @abc.abstractmethod
    def jacobian(self, reference_point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def jacobian_det(self, reference_point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def inverse_jacobian(self, reference_point):
        raise Exception("Abstract method called!")