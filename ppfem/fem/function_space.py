

class FunctionSpace(object):

    def __init__(self, element, domain, mapping_element):
        self._element = element
        self._domain = domain
        # TODO: analyze domain and element an setup function space dofs
        # loop over the corresponding class of entities (edge, face, cell)
        # kann man annehmen, dass alle elemente gleich viele dofs pro knoten zu haben?
        # wie kann man dofs erkennen, die mehreren elementen zugeordnet werden können (hinsichtlich dofs!)
        # '-> support points / dofs aufschreiben?
        # info ob ein dof intern ist oder nicht : ok
        # wenn nicht intern, muss gecheckt werden, ob der physische support point schon aufgenommen wurde

    def function_dim(self):
        pass

    def function_space_dim(self):
        pass
