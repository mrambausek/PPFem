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
import scipy as sp
from ppfem.fem.function import FEFunction

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

    @staticmethod
    @abc.abstractmethod
    def integrate_essential_bc(system_matrix, system_rhs, function_space, indicator_func, bc_func):
        raise Exception("Abstract method called!")


class DefaultSystemAssembler(Assembler):
    def __init__(self):
        Assembler.__init__(self)

    @staticmethod
    def assemble_functionals(discrete_functional, forms, params=None):
        """
        Assembles functionals by adding the local values beginning with the value "discrete_functional"
        :param discrete_functional: a scalar value for the functional to start from
        :param forms: a FormCollection
        :param params: arbitrary params that are handled to the actual local data objects (CellEvalData* and friends)
        :return: the sum over all functionals over their domains (usually a double/float)
        """
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
        """
        Assembles linear forms by adding the local values (arrays) to "discrete_linear_form"
        :param discrete_linear_form: [in/out] an 1d-array-like type that supports element write access via [i]
        :param forms: a FormCollection
        :param params: arbitrary params that are handled to the actual local data objects (CellEvalData* and friends)
        :return: None (The result is written to discrete_linear_form.)
        """
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
        """
        Assembles bilinear forms by adding the local values (matrices) to "discrete_bilinear_form"
        :param discrete_bilinear_form: [in/out] an 2d-array-like type that supports element write access via [i,j]
        :param forms: a FormCollection
        :param params: arbitrary params that are handled to the actual local data objects (CellEvalData* and friends)
        :return: None (The result is written to discrete_bilinear_form.)
        """
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
        """
        Returns the sparsity pattern needed to store all the bilinear_forms contained in forms in a sparse matrix.
        :param forms: a FormCollection
        :return: the sparsity pattern as a 2d-array of which the first row corresponds to row indices and the second
        row to column indices.
        """
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
        return sp.array(global_entries[0], dtype=sp.int64), sp.array(global_entries[1], dtype=sp.int64)

    @staticmethod
    def integrate_essential_bc(system_matrix, system_rhs, function_space, indicator_func, bc_func,
                               current_state=None):
        """
        Integrate simple essential boundary conditions. Simple means bcs like u=u(x).
        Note that the arguments system_matrix, system_rhs and function_space have to be consistent!

        NOTE: This implementation is intented for Bubnov-Galerkin FEM.

        :param system_matrix: [in/out] A (sparse) matrix representing the assembled system matrix.
        :param system_rhs: [in/out] A vector (1d-array) representing the assembled system right-hand-side vector.
        :param function_space: The function space for the constrained quantity.
        :param indicator_func: A function f(x) that returns 'True' for components of the solution is
        constrained at the point 'x'. In case of a scalar solution this function returns a scalar,
        otherwise a list/an array is to be returned,
        :param bc_func: A function u(x) that returns the value of the solution 'u' at the point 'x'.
        In line with indicator_func and the shape of the solution this is either scalar or a list/an array.
        :param current_state: a FEFunction holding the current state of the system. This can be used to
         recompute the essential boundary condition. Especially useful for nonlinear FEM.
        :return: Nothing. System_matrix and system_rhs are manipulated.
        """
        FE_indicator = FEFunction(function_space)
        FE_bc = FEFunction(function_space)

        def _ind(x):
            return indicator_func(x).astype(sp.float64)

        FE_indicator.set_dof_values_from_interpolation(_ind)
        FE_bc.set_dof_values_from_interpolation(bc_func)

        if current_state is not None:
            FE_bc.set_dof_values(FE_bc.dof_values() - current_state.dof_values())

        constrained_indices = sp.where(FE_indicator.dof_values() == 1)
        constrained_dof_values = FE_bc.dof_values()[constrained_indices]
        #print(constrained_dof_values)

        for iv in zip(constrained_indices[0], constrained_dof_values):
            i, value = iv
            for ii in range(system_matrix.shape[0]):
                system_matrix[i, ii] = 0.0
                system_rhs[ii] -= value * system_matrix[ii, i]
                system_matrix[ii, i] = 0.0
            system_matrix[i, i] = 1.0
            system_rhs[i] = value

    @staticmethod
    def _assemble_local_cell_linear_form(local_linear_form, dof_index_array, global_linear_form):
        i = 0
        for I in dof_index_array:
            global_linear_form[I] += local_linear_form[i]
            i += 1

    @staticmethod
    def _assemble_local_cell_bilinear_form(local_bilinear_form, test_dof_index_array, trial_dof_index_array,
                                           global_bilinear_form):
        i = 0
        for I in test_dof_index_array:
            j = 0
            for J in trial_dof_index_array:
                global_bilinear_form[I, J] += local_bilinear_form[i, j]
                j += 1
            i += 1
