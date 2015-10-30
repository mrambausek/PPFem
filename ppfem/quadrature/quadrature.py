import abc


class QPData(object):
    def __init__(self, point, weight):
        self.point = point
        self.weight = weight


class Quadrature(object):

    def __init__(self):
        self._quadrature_data = None

    def __call__(self, quadrature_functor):
        qsum = 0.0
        for qp_data in self._quadrature_data:
            pq = qp_data.point
            wq = qp_data.weight
            qsum += quadrature_functor(pq) * wq

        return qsum


class QuadratureFunctor(abc.ABC):

    def __init__(self, functor, mapping):
        self._functor = functor
        self._mapping = mapping

    def __call__(self, point):
        return self._eval_functor(point) * self._mapping.jacobian_det(point)

    @abc.abstractmethod
    def _eval_functor(self, point):
        raise Exception("Abstract method called!")


class FunctionFunctor(QuadratureFunctor):

    def __init__(self, function):
        self._function = function

    def _eval_functor(self, point):
        real_point = self._mapping.real_point(point)
        return self._function(real_point)