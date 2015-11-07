import abc
import scipy as sp
import scipy.linalg as spl

class Mapping(abc.ABC):

    def __init__(self):
        pass

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


class FEMapping(Mapping):
    def __init__(self, element):
        self._element = element
        self._mesh_entity = None

    def set_mesh_entity(self, mesh_entity):
        self._mesh_entity

    def clear_mesh_entity(self):
        self._mesh_entity = None

    def map_point(self, reference_point):
        return self._element.function_value(self._mesh_entity.vertices(),
                                            reference_point)

    def jacobian(self, reference_point):
        return self._element.function_gradient(self._mesh_entity.vertices(),
                                               reference_point)

    def jacobian_det(self, reference_point):
        J = self.jacobian(reference_point)
        if sp.isscalar(J):
            return 1/sp.array([[J]])
        elif sp.all(sp.array(J.shape) == 1):
            return J.reshape((1,1))
        elif len(J.shape) == 1 or J.shape[0] == 1 or J.shape[1] == 1:
            raise NotImplementedError("Maybe a projection approach is needed here.")
        elif J.shape[0] == J.shape[1]:
            return spl.det(J)
        else:
            raise NotImplementedError("Computing Jacobian \"determinant\" not implemented for shpae=({0:d}, {1:d})"
                                      .format(J.shape))

    def jacobian_inv(self, reference_point):
        J = self.jacobian(reference_point)
        if sp.isscalar(J):
            return sp.array([[1/J]])
        elif sp.all(sp.array(J.shape) == 1):
            return 1/J.reshape((1,1))
        elif len(J.shape) == 1 or J.shape[0] == 1 or J.shape[1] == 1:
            raise NotImplementedError("Maybe a projection approach is needed here.")
        elif J.shape[0] == J.shape[1]:
            return spl.inv(J)
