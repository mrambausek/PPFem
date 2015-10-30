import sympy as sy
from sympy.utilities.autowrap import autowrap


def product(factors, a, b):
    res = sy.sympify(1)
    for i in range(a-1, b):
        res *= factors(i)
    return res


class LagrangeBasis(object):
    def __init__(self, points, index, dimension=1):
        x = sy.symbols('x')
        x_i = points[index][0]

        def f(j):
            return (x - points[j][0]) / (x_i - points[j][0])

        if dimension == 1:
            self._L = product(f, 1, index) * product(f, index+2, len(points))
            self._value = autowrap(self._L)
            self._grad = autowrap(sy.diff(self._L, x), args=(x,))
        else:
            self._L = sy.Matrix([product(f, 1, index) * product(f, index+2, len(points)) for i in range(dimension)])
            self._value = autowrap(self._L, args=[x])
            self._grad = autowrap(
                sy.Matrix([sy.diff(self._L[i], x) for i in range(dimension)]),
                args=[x]
            )

    def L(self):
        return self._L

    def value(self, point):
        return self._value(point[0])

    def gradient(self, point):
        return self._grad(point[0])

