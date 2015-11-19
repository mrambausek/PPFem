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

import abc
import scipy.sparse.linalg as spl


class LinearSolver(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def solve(self, lhs_matrix, lhs_vector, solution_vector):
        raise Exception('Abstract method called!')


class DirectSparseSolver(LinearSolver):

    def __init__(self):
        LinearSolver.__init__(self)

    def solve(self, lhs_matrix, lhs_vector, solution_vector):
        solution_vector[:] = spl.spsolve(lhs_matrix, lhs_vector)
