__author__ = 'mrambausek'


class Function(object):

    def __init__(self, function_space):
        self._function_space = function_space
        self._dof_values = None

    def interpolate_function(self, function):
        raise NotImplementedError("Please implement me.")

    def dof_values(self):
        return self._dof_values

    def set_dof_values(self, new_values):
        # TODO: add some checks
        self._dof_values = new_values
