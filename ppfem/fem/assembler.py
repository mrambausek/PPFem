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
# import scipy as sp
from ppfem.fem.form import LinearForm, BilinearForm


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
    def assemble_linear_forms(self, global_linear_form, global_state, params=None):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_bilinear_forms(self, global_bilinear_form, global_state, params=None):
        raise Exception("Abstract method called!")


class DefaultSystemAssembler(Assembler):
    def __init__(self, forms):
        self._linear_forms = []
        self._bilinear_forms = []
        for f in forms:
            if isinstance(f, LinearForm):
                self._linear_forms.append(f)
            if isinstance(f, BilinearForm):
                self._bilinear_forms.append(f)

    def assemble_linear_forms(self, global_linear_form, global_state, params=None):
        for L in self._linear_forms:
            # FIXME: check what the linear form actually provides!
            # TODO: implement assembly of interior and exterior face terms
            for e in L.test_function_space.mesh_entity_iterator():
                eval_data = L.get_cell_eval_data_linear_form(e, global_state, params)
                local_linear_form = L.local_cell_linear_form(eval_data)
                self._assemble_local_cell_linear_form(local_linear_form, e.index, global_linear_form)

    def assemble_bilinear_forms(self, global_bilinear_form, global_state, params=None):
        for a in self._bilinear_forms:
            # FIXME: check what the bilinear form actually provides!
            # TODO: implement assembly of interior and exterior face terms
            for e in a.trial_function_space.mesh_entity_iterator():
                eval_data = a.get_cell_eval_data_bilinear_form(e, global_state, params)
                local_bilinear_form = a.local_cell_bilinear_form(eval_data)
                self._assemble_local_cell_bilinear_form(local_bilinear_form, e.index, global_bilinear_form)

    def get_sparsity(self):
        global_entries = ([], [])

        for a in self._bilinear_forms:
            for e in a.trial_function_space.mesh_entity_iterator():
                trial_dofs = a.trial_function_space.get_element_dof_index_array(e.index)
                test_dofs = a.test_function_space.get_element_dof_index_array(e.index)
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

    def _assemble_local_cell_linear_form(self, local_linear_form, mesh_entity_index, global_linear_form):
        dofs = self._test_function_space.get_element_dof_index_array(mesh_entity_index)

        i = 0
        for I in dofs:
            global_linear_form[I] = local_linear_form[i]
            i += 1

    def _assemble_local_cell_bilinear_form(self, local_bilinear_form, mesh_entity_index, global_bilinear_form):
        trial_dofs = self._trial_function_space.get_element_dof_index_array(mesh_entity_index)
        test_dofs = self._test_function_space.get_element_dof_index_array(mesh_entity_index)

        i = 0
        for I in test_dofs:
            j = 0
            for J in trial_dofs:
                global_bilinear_form[I, J] = local_bilinear_form[i, j]
                j += 1
            i += 1
