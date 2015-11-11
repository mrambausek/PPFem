import ppfem.quadrature.quadrature as quadrature
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
                    quadrature.QPData(Point( sqrt(1/3), index=1), 1.0)
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

