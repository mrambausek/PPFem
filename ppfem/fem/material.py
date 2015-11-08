import abc


class Material(abc.ABC):

    @abc.abstractmethod
    def compute_rhs(self, physical_coords, function_value, function_gradient):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def compute_lhs(self, physical_coords, function_value, function_gradient):
        raise Exception("Abstract method called!")