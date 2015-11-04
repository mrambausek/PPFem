

class DOF(object):
    def __init__(self, number, value=0.0, vertex=None):
        self.number = number
        self.value = value
        self.associated_vertex = vertex