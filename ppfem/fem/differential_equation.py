import abc


class DifferentialEquation(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def linear_form(self, element, quadrature, params):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def bilinear_form(self, element, quadrature, params):
        raise Exception("Abstract method called!")