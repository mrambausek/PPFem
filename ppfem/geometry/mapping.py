import abc
import numpy as np
import numpy.linalg as npl


class Mapping(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def jacobian_det(self, point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def map_point(self, point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def map_points(self, points):
        return [self.map_point(p) for p in points]

    @abc.abstractmethod
    def jacobian(self, point):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def inverse_jacobian(self, point):
        raise Exception("Abstract method called!")