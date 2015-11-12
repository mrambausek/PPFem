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
import scipy as sp
import scipy.linalg as spl


class Mapping(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def map_point(self, reference_point):
        raise Exception("Abstract method called!")

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
        self._mesh_entity = mesh_entity

    def localize(self, mesh_entity):
        self.set_mesh_entity(mesh_entity)
        return self

    def clear_mesh_entity(self):
        self._mesh_entity = None

    def map_point(self, reference_point):
        return self._element.function_value(self._mesh_entity.vertex_coords(),
                                            reference_point)

    def jacobian(self, reference_point):
        J = self._element.function_gradient(self._mesh_entity.vertex_coords(),
                                            reference_point)
        if sp.isscalar(J):
            return sp.array([[J]])
        elif sp.all(sp.array(J.shape) == 1):
            return J.reshape((1,1))
        else:
            return J

    def jacobian_det(self, reference_point):
        J = self.jacobian(reference_point)
        if sp.isscalar(J):
            return J
        elif sp.all(sp.array(J.shape) == 1):
            # an array with only one entry
            return sp.asscalar(J)
        elif len(J.shape) == 1 or J.shape[0] == 1 or J.shape[1] == 1:
            return sp.linalg.norm(J)
        elif J.shape[0] == J.shape[1]:
            # a square matrix
            return spl.det(J)
        else:
            raise NotImplementedError("Computing Jacobian \"determinant\" not implemented for shpae=({0:d}, {1:d})"
                                      .format(J.shape))

    def inverse_jacobian(self, reference_point):
        J = self.jacobian(reference_point)
        if sp.isscalar(J):
            return sp.array([[1/J]])
        elif sp.all(sp.array(J.shape) == 1):
            return 1/J.reshape((1,1))
        elif len(J.shape) == 1 or J.shape[0] == 1 or J.shape[1] == 1:
            raise NotImplementedError("Maybe a projection approach is needed here.")
        elif J.shape[0] == J.shape[1]:
            print(J.shape, J)
            return spl.inv(J)
