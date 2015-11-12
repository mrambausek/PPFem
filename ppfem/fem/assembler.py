# PPFem: A educational finite element code
# Copyright (C) 2015  Matthias Rambausek
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

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
    def get_sparsity(self):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_rhs(self, global_rhs, params):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def assemble_lhs(self, global_lhs, params):
        raise Exception("Abstract method called!")
