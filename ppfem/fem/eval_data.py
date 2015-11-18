# PPFem: An educational finite element code
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


class EvalData(object):
    def __init__(self, attr_dict):
        self.add_attributes(attr_dict)

    def add_attributes(self, attr_dict):
        for k in attr_dict.keys():
            self.add_attribute(k, attr_dict[k])

    def add_attribute(self, name, value):
        self.__setattr__(name, value)

    def remove_attribute(self, name):
        self.__delattr__(name)

    def remove_attributes(self, names):
        for name in names:
            self.remove_attribute(name)
