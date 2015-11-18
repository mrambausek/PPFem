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

    @staticmethod
    @abc.abstractmethod
    def get_sparsity(self, forms):
        raise Exception("Abstract method called!")

    @staticmethod
    @abc.abstractmethod
    def assemble_functionals(discrete_functional, forms, params=None):
        raise Exception("Abstract method called!")

    @staticmethod
    @abc.abstractmethod
    def assemble_linear_forms(discrete_linear_form, forms, params=None):
        raise Exception("Abstract method called!")

    @staticmethod
    @abc.abstractmethod
    def assemble_bilinear_forms(discrete_bilinear_form, forms, params=None):
        raise Exception("Abstract method called!")


class DefaultSystemAssembler(Assembler):
    def __init__(self):
        Assembler.__init__(self)

    @staticmethod
    def assemble_functionals(discrete_functional, forms, params=None):
        for f in forms.functional_iterator():
            # FIXME: check what the linear form actually provides!
            # TODO: implement assembly of interior and exterior face terms
            for e in f.mesh_entity_iterator():
                eval_data = f.get_cell_eval_data_functional(e, params)
                local_functional = f.local_cell_functional(eval_data)
                discrete_functional += local_functional
        return discrete_functional

    @staticmethod
    def assemble_linear_forms(discrete_linear_form, forms, params=None):
        for L in forms.linear_form_iterator():
            # FIXME: check what the linear form actually provides!
            # TODO: implement assembly of interior and exterior face terms
            for e in L.mesh_entity_iterator():
                eval_data = L.get_cell_eval_data_linear_form(e, params)
                local_linear_form = L.local_cell_linear_form(eval_data)
                DefaultSystemAssembler._assemble_local_cell_linear_form(
                    local_linear_form,
                    L.test_function_space.get_element_dof_index_array(e.index),
                    discrete_linear_form
                )

    @staticmethod
    def assemble_bilinear_forms(discrete_bilinear_form, forms, params=None):
        for a in forms.bilinear_form_iterator():
            # FIXME: check what the bilinear form actually provides!
            # TODO: implement assembly of interior and exterior face terms
            for e in a.mesh_entity_iterator():
                eval_data = a.get_cell_eval_data_bilinear_form(e, params)
                local_bilinear_form = a.local_cell_bilinear_form(eval_data)
                DefaultSystemAssembler._assemble_local_cell_bilinear_form(
                    local_bilinear_form,
                    a.test_function_space.get_element_dof_index_array(e.index),
                    a.trial_function_space.get_element_dof_index_array(e.index),
                    discrete_bilinear_form
                )

    def get_sparsity(self, forms):
        global_entries = ([], [])

        for a in forms.bilinear_form_iterator():
            for e in a.mesh_entity_iterator():
                trial_dofs = a.trial_function_space.get_element_dof_index_array(e.index)
                test_dofs = a.test_function_space.get_element_dof_index_array(e.index)
                i = 0
                for I in test_dofs:
                    j = 0
                    for J in trial_dofs:
                        p = [I, J]
                        if not find_index_pair(p, global_entries):
                            global_entries[0].append(I)
                            global_entries[1].append(J)
                        j += 1
                    i += 1
        return global_entries

    @staticmethod
    def _assemble_local_cell_linear_form(local_linear_form, dof_index_array, global_linear_form):
        i = 0
        for I in dof_index_array:
            global_linear_form[I] = local_linear_form[i]
            i += 1

    @staticmethod
    def _assemble_local_cell_bilinear_form(local_bilinear_form, test_dof_index_array, trial_dof_index_array,
                                           global_bilinear_form):
        i = 0
        for I in test_dof_index_array:
            j = 0
            for J in trial_dof_index_array:
                global_bilinear_form[I, J] = local_bilinear_form[i, j]
                j += 1
            i += 1
