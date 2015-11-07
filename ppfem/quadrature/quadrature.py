import abc


class QPData(object):
    def __init__(self, point, weight):
        self.point = point
        self.weight = weight


class Quadrature(abc.ABC):

    def __init__(self):
        self._quadrature_data = None

    def number_of_quadrature_points(self):
        return len(self._quadrature_data)

    def quadrature_data(self):
        return self._quadrature_data

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
        return self._functor(point) * self._mapping.jacobian_det(point)
