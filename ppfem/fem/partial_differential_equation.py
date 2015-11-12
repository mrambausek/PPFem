import abc
from ppfem.fem.assembler import find_index_pair, Assembler


class PDE(Assembler):
    def __init__(self, test_function_space, trial_function_space, quadrature, fe_functions=None):
        self._test_function_space = test_function_space
        self._trial_function_space = trial_function_space
        self._quadrature = quadrature
        self.number_of_trial_dofs = trial_function_space.number_of_dofs
        self.number_of_test_dofs = test_function_space.number_of_dofs
        self._local_test_function_space = None
        self._local_trial_function_space = None

        if fe_functions is None:
            self._fe_functions = []
        else:
            self._fe_functions = fe_functions

    def localize(self, mesh_entity):
        self._local_test_function_space = self._test_function_space.localize(mesh_entity)
        self._local_trial_function_space = self._trial_function_space.localize(mesh_entity)
        for f in self._fe_functions:
            f.localize(mesh_entity)

    def assemble_rhs(self, global_rhs, params):
        for e in self._trial_function_space.mesh_entity_iterator():
            self.localize(e)
            f = self.local_linear_form(params)
            self._assemble_local_linear_form(f, e.index, global_rhs)

    def assemble_lhs(self, global_lhs, params):
        for e in self._trial_function_space.mesh_entity_iterator():
            self.localize(e)
            a = self.local_bilinear_form(params)
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

    @abc.abstractmethod
    def local_linear_form(self, params):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_bilinear_form(self, params):
        raise Exception("Abstract method called!")
