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

import ppfem.fem.form as form


class PDE(form.LinearForm, form.BilinearForm):
    """A class implementing LinearForm and BilinearForm."""
    def __init__(self, test_function_space, trial_function_space, quadrature, fe_functions=None, non_fe_functions=None,
                 mapping=None, mesh=None, subdomain=None):

        form.LinearForm.__init__(self, test_function_space, quadrature, fe_functions=fe_functions,
                                 non_fe_functions=non_fe_functions, mapping=mapping, mesh=mesh, subdomain=subdomain)

        form.BilinearForm.__init__(self, test_function_space, trial_function_space, quadrature,
                                   fe_functions=fe_functions, non_fe_functions=non_fe_functions, mapping=mapping,
                                   mesh=mesh, subdomain=subdomain)
