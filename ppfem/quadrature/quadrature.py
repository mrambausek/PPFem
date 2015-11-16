# PPFem: A educational finite element code
# Copyright (C) 2015  Matthias Rambausek
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import abc


class QPData(object):
    def __init__(self, point, weight):
        self.point = point
        self.weight = weight


class Quadrature(abc.ABC):

    def __init__(self):
        self._quadrature_data = []

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

    def __call__(self, q_point):
        return self._functor(q_point) * self._mapping.jacobian_det(q_point)
