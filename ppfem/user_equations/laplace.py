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

from ppfem.fem.partial_differential_equation import PDE
import scipy as sp


class Laplace(PDE):

    def __init__(self, test_function_space, trial_function_space, quadrature, dim=1):
        PDE.__init__(self, test_function_space, trial_function_space, quadrature)
        self._dim = dim

    def local_linear_form(self, params):
        raise NotImplementedError("Implement me!")

    def local_bilinear_form(self, params):
        raise NotImplementedError("Implement me!")
