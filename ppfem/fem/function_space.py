

class FunctionSpace(object):

    def __init__(self, mesh, element_map):
        self._element_map = element_map
        self._mesh = mesh

    def function_dim(self):
        pass

    def function_space_dim(self):
        pass
