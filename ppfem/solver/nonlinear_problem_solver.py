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

from ppfem.fem.assembler import get_essential_bc_data
from ppfem.solver import linear_solver as linsolve


class SimpleNewtonSolver(object):

    def __init__(self, forms, assembler, state, bcs=[], linear_solver=linsolve.DirectSparseSolver()):
        """
        A very simple, but for some linear(ized) problems suffient solver routine.
        The constructor already does some sanity checks and sets up the storage (matrices etc.)
        :param forms: a FormCollection containing linear and bilinear forms.
        :param assembler: an Assembler
        :param state: an FEFunction representing the current state.
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

        if self._V_trial != state.function_space:
            raise Exception("The function space of the state must be the same as for the trial function.")

        sparsity = self._assembler.get_sparsity(self._forms)
        initial_data = sp.zeros_like(sparsity[0], dtype=sp.float64)
        self.lhs_matrix = sparse.csr_matrix((initial_data, sparsity))
        self.rhs_vector = sp.zeros(self.lhs_matrix.shape[1])
        self.solution = sp.zeros_like(self.rhs_vector)
        self.state = state
        self.info = {'converged': False, 'iter': 0}
        self._state_backup = None
        self._bc_indices = None
        self._free_indices = None
        self._bc_values = None

    def _backup_state(self):
        self._state_backup = self.state.dof_values()

    def _restore_state(self):
        self.state.set_dof_values(self._state_backup)

    def _test_convergence(self, rhs_tol, dx_tol, bc_tol):
        self.info['||rhs||'] = sp.linalg.norm(self.rhs_vector[self._free_indices])
        self.info['||dx||'] = sp.linalg.norm(self.solution)
        self.info['||bc||'] = sp.linalg.norm(self._bc_values)
        if self.info['||rhs||'] <= rhs_tol \
            and self.info['||dx||'] <= dx_tol \
            and self.info['||bc||'] <= bc_tol:
                self.info['converged'] = True
        return self.info['converged']

    def print_header(self):
        print("{:>4s} -- {:^9s} -- {:^9s} -- {:^9s} -- {:^3s}".format('#It:','||RHS||','||dx||','||dx_bc||', 'CON'))

    def print_info(self):
        print("{niter:>4d} -- {rhs:7.3e} -- {dx:7.3e} -- {bc:7.3e} -- {converged:^3b}".format(
                                                                         niter=self.info['iter'],
                                                                         rhs=self.info['||rhs||'],
                                                                         dx=self.info['||dx||'],
                                                                         bc=self.info['||bc||'],
                                                                         converged=self.info['converged'])
        )

    def _bc_data(self):
        #TODO: many optimizations possible here!
        indices = sp.empty(0, dtype=sp.int64)
        values = sp.empty(0)
        for bc in self._bcs:
            ind, val = get_essential_bc_data(bc['indicator'], bc['func'], self.state.function_space, self.state)
            indices = sp.r_[indices,ind[0]]
            values = sp.r_[values, val]
        self._bc_indices = indices
        self._free_indices = sp.array(list(filter(lambda i: i not in indices,
                                             range(self.state.function_space.number_of_dofs))))
        self._bc_values = values

    def _assemble_linear_forms(self, params=None):
        self.rhs_vector[:] = 0.0
        self._assembler.assemble_linear_forms(self.rhs_vector, self._forms, params=params)

    def _assemble_bilinear_forms(self, params=None):
        self.lhs_matrix.data[:] = 0.0
        self._assembler.assemble_bilinear_forms(self.lhs_matrix, self._forms, params=params)

    def solve(self, rhs_tol=1e-8, bc_tol=1e-8, dx_tol=1e-8, max_iter=10, params=None):
        """
        The actual solving of the linear system happens here.
        :param dx_tol: tolerance for ||dx|| (convergence criterion)
        :param rhs_tol: tolerance for ||RHS|| (convergence criterion)
        :param max_iter: max number of iterations
        :param params: params to be handed over to the assembler routine.
        :return:the solution (an array) or the updated state (an FEFunction) if given as keyword argument.
        """
        state = self.state
        self._bc_data()

        self._assemble_linear_forms(params=params)

        self.print_header()

        if self._test_convergence(rhs_tol, dx_tol, bc_tol):
            print("Initial residual and constraint violation already pass the convergence test!")
            return state

        for i in range(max_iter):
            self._assemble_bilinear_forms(params=params)

            for bc in self._bcs:
                self._assembler.integrate_essential_bc(self.lhs_matrix, self.rhs_vector, self._V_trial,
                                                       bc['indicator'], bc['func'],
                                                       current_state=state)

            self._linear_solver.solve(self.lhs_matrix, self.rhs_vector, self.solution)
            self.info['iter'] = i+1

#            self._backup_state()
            state.set_dof_values(state.dof_values() + self.solution)
            self._bc_data()
            self._assemble_linear_forms(params=params)

            if self._test_convergence(rhs_tol, dx_tol, bc_tol):
                self.print_info()
                break
            self.print_info()

        return state
