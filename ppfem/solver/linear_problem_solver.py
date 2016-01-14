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

import scipy as sp
import scipy.sparse as sparse
from ppfem.fem import linear_solver as linsolve


class SimpleLinearProblemSolver(object):

    def __init__(self, forms, assembler, bcs=[], linear_solver=linsolve.DirectSparseSolver()):
        """
        A very simple, but for some linear(ized) problems suffient solver routine.
        The constructor already does some sanity checks and sets up the storage (matrices etc.)
        :param forms: a FormCollection containing linear and bilinear forms.
        :param assembler: an Assembler
        :param bcs: an iterable of dicts like {'indicator': bc_indicator, 'func': bc_func}.
        Have a look at the docstrings of "integrate_essential_bc"-routines in assembler.py.
        :param linear_solver: a LinearSolver for sparse matrices
        """
        self._forms = forms
        self._assembler = assembler
        self._linear_solver = linear_solver
        self._bcs = bcs

        _test = self._forms.test_function_spaces()
        _trial = self._forms.trial_function_spaces()

        if len(_test) > 1:
            raise Exception("All forms must use the same test function spaces.")
        if len(_trial) > 1:
            raise Exception("All forms must use the same trial function spaces.")

        self._V_test = _test[0]
        self._V_trial = _trial[0]

        sparsity = self._assembler.get_sparsity(self._forms)
        initial_data = sp.zeros_like(sparsity[0], dtype=sp.float64)
        self._lhs_matrix = sparse.csr_matrix((initial_data, sparsity))
        self._rhs_vector = sp.zeros(self._lhs_matrix.shape[1])
        self._solution = sp.zeros_like(self._rhs_vector)

    def solve(self, state=None, params=None):
        """
        The actual solving of the linear system happens here.
        :param params: params to be handed over to the assembler routine.
        :param state: an FEFunction representing the current state.
        :return:the solution (an array) or the updated state (an FEFunction) if given as keyword argument.
        """
        self._lhs_matrix.data[:] = 0.0
        self._rhs_vector[:] = 0.0

        self._assembler.assemble_bilinear_forms(self._lhs_matrix, self._forms, params=params)
        self._assembler.assemble_linear_forms(self._rhs_vector, self._forms, params=params)
        for bc in self._bcs:
            self._assembler.integrate_essential_bc(self._lhs_matrix, self._rhs_vector, self._V_trial,
                                                   bc['indicator'], bc['func'],
                                                   current_state=state)

        self._linear_solver.solve(self._lhs_matrix, self._rhs_vector, self._solution)
        # print()
        # print(self._lhs_matrix.toarray())
        # print(self._rhs_vector)
        # print(self._solution)
        # print()
        if state is not None:
            state.set_dof_values(state.dof_values() + self._solution)
            return state
        else:
            return self._solution
