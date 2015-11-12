import abc


def find_index_pair(pair, data):
    d_i = data[0]
    d_j = data[1]
    for i in range(len(data[0])):
        p_global = [d_i[i], d_j[i]]

        if p_global[0] == pair[0] and p_global[1] == pair[1]:
            return True
    return False


class Assembler(abc.ABC):
    def __init__(self):
        pass

    @abc.abstractmethod
    def get_sparsity(self, sp_format='coo', mesh_entities=None):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_rhs(self, global_rhs):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_lhs(self, global_lhs):
        raise Exception("Abstract method called!")
