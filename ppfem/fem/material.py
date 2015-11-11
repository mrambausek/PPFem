import abc


class ContractionData(object):
    func = 'func'
    grad = 'grad'

    def __init__(self, rhs_contraction, lhs_contraction):
        self.rhs = rhs_contraction
        self.lhs = lhs_contraction


class Material(abc.ABC):
    def __init__(self):
        pass

    @staticmethod
    @abc.abstractmethod
    def contraction_data():
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def compute_rhs(self, physical_coords, function_value, function_gradient):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def compute_lhs(self, physical_coords, function_value, function_gradient):
        raise Exception("Abstract method called!")