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


class CellEvalDataBase(object):
    def __init__(self, local_state, mapping, quadrature, params):
        self.local_state = local_state
        self.mapping = mapping
        self.quadrature = quadrature
        self.params = params


class CellEvalDataLinearForm(CellEvalDataBase):
    def __init__(self, local_test_space, local_state, mapping, quadrature, params):
        CellEvalDataBase.__init__(self, local_state, mapping, quadrature, params)
        self.local_test_space = local_test_space


class CellEvalDataBilinearForm(CellEvalDataLinearForm):
    def __init__(self, local_test_space, local_trial_space, local_state, mapping, quadrature, params):
        CellEvalDataLinearForm.__init__(self, local_test_space, local_state, mapping, quadrature, params)
        self.local_trial_space = local_trial_space


# FIXME: *FaceEvalData* not tested and not a stable design yet!
class InteriorFaceEvalDataBase(object):
    def __init__(self, local_states, mappings, quadrature, params):
        self.local_states = local_states
        self.mapping = mappings
        self.quadrature = quadrature
        self.params = params


class InteriorFaceEvalDataLinearForm(InteriorFaceEvalDataBase):
    def __init__(self, local_test_spaces, local_states, mappings, quadrature, params):
        InteriorFaceEvalDataBase.__init__(self, local_states, mappings, quadrature, params)
        self.local_test_spaces = local_test_spaces


class InteriorFaceEvalDataBilinearForm(InteriorFaceEvalDataLinearForm):
    def __init__(self, local_test_spaces, local_trial_spaces, local_states, mappings, quadrature, params):
        InteriorFaceEvalDataLinearForm.__init__(self, local_test_spaces, local_states, mappings, quadrature, params)
        self.local_trial_spaces = local_trial_spaces


class ExteriorFaceEvalDataBase(object):
    def __init__(self, local_state, mapping, quadrature, params):
        self.local_state = local_state
        self.mapping = mapping
        self.quadrature = quadrature
        self.params = params


class ExteriorFaceEvalDataLinearForm(ExteriorFaceEvalDataBase):
    def __init__(self, local_test_space, local_state, mapping, quadrature, params):
        ExteriorFaceEvalDataBase.__init__(self, local_state, mapping, quadrature, params)
        self.local_test_space = local_test_space


class ExteriorFaceEvalDataBilinearForm(ExteriorFaceEvalDataLinearForm):
    def __init__(self, local_test_space, local_trial_space, local_state, mapping, quadrature, params):
        ExteriorFaceEvalDataLinearForm.__init__(self, local_test_space, local_state, mapping, quadrature, params)
        self.local_trial_space = local_trial_space


class Form(abc.ABC):
    cells = "cells"
    interior_faces = "interior_faces"
    exterior_faces = "exterior_faces"

    def __init__(self, quadrature):
        self.quadrature = quadrature

    @abc.abstractmethod
    def implements_quadrature_on(self, entity_type=None):
        """
        Implementation should give information whether a cell_linear_form etc is provided or not.
        :param entity_type: Form.cell, .interior_faces, .exterior_faces to query if such integrals are computed.
        """
        raise Exception("Abstract method called!")


class LinearForm(Form):
    def __init__(self, test_function_space, quadrature):
        Form.__init__(self, quadrature)
        self.test_function_space = test_function_space
        self.mapping = self.test_function_space.get_mapping()

    @abc.abstractmethod
    def local_cell_linear_form(self, cell_eval_data_linear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_interior_face_linear_form(self, interior_face_eval_data_linear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_exterior_face_linear_form(self, exterior_face_eval_data_linear_form):
        raise Exception("Abstract method called!")

    def get_cell_eval_data_linear_form(self, mesh_entity, solution, params=None, mapping=None):
        if mapping is None:
            mapping = self.mapping
        return CellEvalDataLinearForm(self.test_function_space.localize(mesh_entity),
                                      solution.localize(mesh_entity),
                                      self.quadrature,
                                      mapping.localize(mesh_entity),
                                      params)

    def get_interior_face_eval_data_linear_form(self, mesh_entities, solution, params=None, mapping=None):
        if mapping is None:
            mapping = self.mapping
        return InteriorFaceEvalDataLinearForm(
            [self.test_function_space.localize(mesh_entities[0]),
             self.test_function_space.localize(mesh_entities[1])],
            [solution.localize(mesh_entities[0]),
             solution.localize(mesh_entities[1])],
            self.quadrature,
            [mapping.localize(mesh_entities[0]),
             mapping.localize(mesh_entities[1])],
            params)

    def get_exterior_face_eval_data_linear_form(self, mesh_entity, solution, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        if mapping is None:
            mapping = self.mapping
        return ExteriorFaceEvalDataLinearForm(self.test_function_space.localize(mesh_entity),
                                              solution.localize(mesh_entity),
                                              self.quadrature,
                                              mapping.localize(mesh_entity),
                                              params)

    def get_local_size(self, mesh_entity):
        return self.test_function_space.get_number_of_global_element_dofs(mesh_entity)

    def get_test_function_space_dim(self):
        return self.test_function_space.function_space_dim()


class BilinearForm(Form):
    def __init__(self, test_function_space, trial_function_space, quadrature):
        Form.__init__(self, quadrature)
        self.test_function_space = test_function_space
        self.trial_function_space = trial_function_space
        self.mapping = self.test_function_space.get_mapping()

    @abc.abstractmethod
    def local_cell_bilinear_form(self, cell_eval_data_bilinear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_interior_face_bilinear_form(self, interior_face_eval_data_bilinear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_exterior_face_bilinear_form(self, exterior_face_eval_data_bilinear_form):
        raise Exception("Abstract method called!")

    def get_cell_eval_data_bilinear_form(self, mesh_entity, solution, params=None, mapping=None):
        if mapping is None:
            mapping = self.mapping
        return CellEvalDataBilinearForm(self.test_function_space.localize(mesh_entity),
                                        self.trial_function_space.localize(mesh_entity),
                                        solution.localize(mesh_entity),
                                        self.quadrature,
                                        mapping.localize(mesh_entity),
                                        params)

    def get_interior_face_eval_data_bilinear_form(self, mesh_entities, solution, params=None, mapping=None):
        if mapping is None:
            mapping = self.mapping
        return InteriorFaceEvalDataBilinearForm(
            [self.test_function_space.localize(mesh_entities[0]),
             self.test_function_space.localize(mesh_entities[1])],
            [self.trial_function_space.localize(mesh_entities[0]),
             self.trial_function_space.localize(mesh_entities[1])]
            [solution.localize(mesh_entities[0]),
             solution.localize(mesh_entities[1])],
            self.quadrature,
            [mapping.localize(mesh_entities[0]),
             mapping.localize(mesh_entities[1])],
            params)

    def get_exterior_face_eval_data_bilinear_form(self, mesh_entity, solution, params=None, mapping=None):
        if mapping is None:
            mapping = self.mapping
        return ExteriorFaceEvalDataBilinearForm(self.test_function_space.localize(mesh_entity),
                                                self.trial_function_space.localize(mesh_entity),
                                                solution.localize(mesh_entity),
                                                self.quadrature,
                                                mapping.localize(mesh_entity),
                                                params)

    def get_local_shape(self, mesh_entity):
        return (self.test_function_space.get_number_of_global_element_dofs(mesh_entity),
                self.trial_function_space.get_number_of_global_element_dofs(mesh_entity))

    def get_test_function_space_dim(self):
        return self.test_function_space.function_space_dim()

    def get_trial_function_space_dim(self):
        return self.trial_function_space.function_space_dim()
