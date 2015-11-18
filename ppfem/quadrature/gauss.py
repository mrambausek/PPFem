# PPFem: An educational finite element code
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

from . import quadrature
from ppfem.geometry.point import Point
from math import sqrt


class QGauss(quadrature.Quadrature):

    def __init__(self, shape, degree):
        quadrature.Quadrature.__init__(self)
        if shape == "line":
            if degree <= 1:
                self._quadrature_data = [
                    quadrature.QPData(Point(0.0, index=0), 2.0)
                ]
            elif degree <= 3:
                self._quadrature_data = [
                    quadrature.QPData(Point(-sqrt(1/3), index=0), 1.0),
                    quadrature.QPData(Point(sqrt(1/3), index=1), 1.0)
                ]
            elif degree <= 5:
                self._quadrature_data = [
                    quadrature.QPData(Point(-sqrt(3/5), index=0),
                                      5/9),
                    quadrature.QPData(Point(0.0, index=1),
                                      8/9),
                    quadrature.QPData(Point(sqrt(3/5), index=2),
                                      5/9)
                ]
            else:
                raise Exception("Gauss quadrature not implemented for degree \"{0:d}\"".format(degree))
        else:
            raise Exception("Gauss quadrature not implemented for shape \"{0:s}\"".format(shape))
