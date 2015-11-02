import abc


class NonlinearProblemSolver(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def solve(self, nonlinear_problem):
        raise Exception('abstract method called!')