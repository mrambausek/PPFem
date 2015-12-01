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

import ppfem.user_elements
import ppfem.user_equations
import ppfem.quadrature

from ppfem.user_elements import *
from ppfem.user_equations import *
from ppfem.quadrature import *

from ppfem.mesh.mesh import Mesh
from ppfem.geometry import Point, Vertex, Line, Face, Cell, Mapping
from ppfem.fem.assembler import DefaultSystemAssembler
from ppfem.fem.form import Functional, LinearForm, BilinearForm, FormCollection
from ppfem.fem.function import FEFunction, FunctionEvaluator
from ppfem.fem.function_space import FunctionSpace
from ppfem.fem.partial_differential_equation import PDE

__all__ = ["Mesh", "Point", "Line", "Vertex", "Face", "Cell", "Mapping", "FunctionSpace", "Functional",
           "LinearForm", "BilinearForm", "FormCollection", "DefaultSystemAssembler", "FEFunction", "FunctionEvaluator",
           "PDE"]

__all__ += ppfem.user_elements.__all__ + ppfem.quadrature.__all__ + ppfem.user_equations.__all__
