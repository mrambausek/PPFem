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


def find_index_pair(pair, data):
    d_i = data[0]
    d_j = data[1]
    for i in range(len(data[0])):
        p_global = [d_i[i], d_j[i]]

        if p_global[0] == pair[0] and p_global[1] == pair[1]:
            return True
    return False


class Assembler(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def get_sparsity(self):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_rhs(self, global_rhs, global_state, params):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_lhs(self, global_lhs, global_state, params):
        raise Exception("Abstract method called!")


class EvalData(object):
    def __init__(self, local_test_space, local_trial_space, local_state, mapping, quadrature):
        self.local_test_space = local_test_space
        self.local_trial_space = local_trial_space
        self.local_state = local_state
        self.mapping = mapping
        self.quadrature = quadrature


class DefaultAssembler(Assembler):
    def __init__(self, test_function_space, trial_function_space, quadrature, mapping=None):
        self._test_function_space = test_function_space
        self._trial_function_space = trial_function_space
        self._quadrature = quadrature
        self.number_of_trial_dofs = trial_function_space.number_of_dofs
        self.number_of_test_dofs = test_function_space.number_of_dofs
        if mapping is not None:
            self._mapping = mapping
        else:
            self._mapping = self._trial_function_space.get_mapping()

    def get_eval_data(self, solution, mesh_entity):
        return EvalData(self._test_function_space.localize(mesh_entity),
                        self._trial_function_space.localize(mesh_entity),
                        solution.localize(mesh_entity),
                        self._quadrature,
                        self._mapping.localize(mesh_entity))

    def assemble_rhs(self, global_rhs, global_state, params):
        for e in self._trial_function_space.mesh_entity_iterator():
            eval_data = self.get_eval_data(global_state, e)
            f = self.local_linear_form(eval_data, params)
            self._assemble_local_linear_form(f, e.index, global_rhs)

    def assemble_lhs(self, global_lhs, global_state, params):
        for e in self._trial_function_space.mesh_entity_iterator():
            eval_data = self.get_eval_data(global_state, e)
            a = self.local_bilinear_form(eval_data, params)
            self._assemble_local_bilinear_form(a, e.index, global_lhs)

    def get_sparsity(self):
        global_entries = ([], [])

        for e in self._trial_function_space.mesh_entity_iterator():
            trial_dofs = self._trial_function_space.get_element_dof_index_array(e.index)
            test_dofs = self._test_function_space.get_element_dof_index_array(e.index)
            i = 0
            for I in test_dofs:
                j = 0
                for J in trial_dofs:
                    p = [I,J]
                    if not find_index_pair(p, global_entries):
                        global_entries[0].append(I)
                        global_entries[1].append(J)
                    j += 1
                i += 1
        return global_entries

    def _assemble_local_linear_form(self, f, mesh_entity_index, global_rhs):
        dofs = self._test_function_space.get_element_dof_index_array(mesh_entity_index)

        i = 0
        for I in dofs:
            global_rhs[I] = f[i]
            i += 1

    def _assemble_local_bilinear_form(self, a, mesh_entity_index, global_lhs):
        trial_dofs = self._trial_function_space.get_element_dof_index_array(mesh_entity_index)
        test_dofs = self._test_function_space.get_element_dof_index_array(mesh_entity_index)

        i = 0
        for I in test_dofs:
            j = 0
            for J in trial_dofs:
                global_lhs[I, J] = a[i, j]
                j += 1
            i += 1
