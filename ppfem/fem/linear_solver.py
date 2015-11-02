import abc
import scipy.sparse as sparse
import scipy.sparse.linalg as spl


class LinearSolver(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def solve(self, lhs_matrix, lhs_vector, solution_vector):
        raise Exception('Abstract method called!')


class DirectSparseSolver(LinearSolver):

    def __init__(self):
        pass

    def solve(self, lhs_matrix, lhs_vector, solution_vector):
        solution_vector = spl.spsolve(lhs_matrix, lhs_vector)
