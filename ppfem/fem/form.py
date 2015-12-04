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

import abc


class CellEvalDataBase(object):
    """
    A simple struct-like type aggregating objects/data for integration of Forms on element-level.
    For the integration of *linear forms, have a look at the appropriate classes below.

    The intention is that instances of this class hold purely local (element-level) data.
    Thus, an assembler may ask a form to provide such data for that it can be used to compute
    the local expression.

    From a design perspective this functionaly could also be totally hidden from the assembler. This might change
    in the future.
    """
    def __init__(self, local_fe_functions, local_non_fe_functions, mapping, quadrature, params):
        self.local_fe_functions = local_fe_functions
        self.local_non_fe_functions = local_non_fe_functions
        self.mapping = mapping
        self.quadrature = quadrature
        self.params = params


class CellEvalDataLinearForm(CellEvalDataBase):
    """
    This class extends its parent to hold a localized function space (local_test_space).
    """
    def __init__(self, local_test_space, local_fe_functions, local_non_fe_functions, mapping, quadrature, params):
        CellEvalDataBase.__init__(self, local_fe_functions, local_non_fe_functions, mapping, quadrature, params)
        self.local_test_space = local_test_space


class CellEvalDataBilinearForm(CellEvalDataLinearForm):
    """
    This class extends its parent to hold one more localized function space (local_trial_space).
    """
    def __init__(self, local_test_space, local_trial_space, local_fe_functions, local_non_fe_functions, mapping,
                 quadrature, params):
        CellEvalDataLinearForm.__init__(self, local_test_space, local_fe_functions, local_non_fe_functions, mapping,
                                        quadrature, params)
        self.local_trial_space = local_trial_space


# FIXME: *FaceEvalData* not tested and not a stable design yet!
class InteriorFaceEvalDataBase(object):
    def __init__(self, local_fe_functions, local_non_fe_functions, mappings, quadrature, params):
        self.local_fe_functions = local_fe_functions
        self.local_non_fe_functions = local_non_fe_functions
        self.mapping = mappings
        self.quadrature = quadrature
        self.params = params


class InteriorFaceEvalDataLinearForm(InteriorFaceEvalDataBase):
    def __init__(self, local_test_spaces, local_fe_functions, local_non_fe_functions, mappings, quadrature, params):
        InteriorFaceEvalDataBase.__init__(self, local_fe_functions, local_non_fe_functions, mappings, quadrature,
                                          params)
        self.local_test_spaces = local_test_spaces


class InteriorFaceEvalDataBilinearForm(InteriorFaceEvalDataLinearForm):
    def __init__(self, local_test_spaces, local_trial_spaces, local_fe_functions, local_non_fe_functions, mappings,
                 quadrature, params):
        InteriorFaceEvalDataLinearForm.__init__(self, local_test_spaces, local_fe_functions, local_non_fe_functions,
                                                mappings, quadrature, params)
        self.local_trial_spaces = local_trial_spaces


class ExteriorFaceEvalDataBase(object):
    def __init__(self, local_fe_functions, local_non_fe_functions, mapping, quadrature, params):
        self.local_fe_functions = local_fe_functions
        self.local_non_fe_functions = local_non_fe_functions
        self.mapping = mapping
        self.quadrature = quadrature
        self.params = params


class ExteriorFaceEvalDataLinearForm(ExteriorFaceEvalDataBase):
    def __init__(self, local_test_space, local_fe_functions, local_non_fe_functions, mapping, quadrature, params):
        ExteriorFaceEvalDataBase.__init__(self, local_fe_functions, local_non_fe_functions, mapping, quadrature, params)
        self.local_test_space = local_test_space


class ExteriorFaceEvalDataBilinearForm(ExteriorFaceEvalDataLinearForm):
    def __init__(self, local_test_space, local_trial_space, local_fe_functions, local_non_fe_functions, mapping,
                 quadrature, params):
        ExteriorFaceEvalDataLinearForm.__init__(self, local_test_space, local_fe_functions, local_non_fe_functions,
                                                mapping, quadrature, params)
        self.local_trial_space = local_trial_space


class Form(abc.ABC):
    """
    An abstract class being the base for implementations of functional, linear forms and bilinear forms.
    This base class provides some very basic interface that is extended by Functional, LinearForm and BilinearForm.
    """

    # some "constants"
    cells = "cells"
    interior_faces = "interior_faces"
    exterior_faces = "exterior_faces"

    def __init__(self, quadrature, fe_functions=None, non_fe_functions=None, mapping=None, mesh=None, subdomain=None):
        self.quadrature = quadrature
        self._mesh = mesh
        self._subdomain = subdomain

        if fe_functions is None:
            self.fe_functions = {}
        else:
            self.fe_functions = fe_functions
        self._localized_fe_functions = {}

        if non_fe_functions is None:
            self.non_fe_functions = {}
        else:
            self.non_fe_functions = non_fe_functions
        self._localized_non_fe_functions = {}

        self.mapping = mapping

        if self.mapping is None:
            for fef in self.fe_functions.values():
                if hasattr(fef, "get_mapping"):
                    try:
                        self.mapping = fef.get_mapping()
                        break
                    except TypeError:
                        pass

        if self._mesh is not None:
            self.set_mesh(mesh, subdomain)

    def _localize_fe_functions(self, mesh_entity):
        for k in self.fe_functions.keys():
            self._localized_fe_functions[k] = self.fe_functions[k].localize(mesh_entity)
        return self._localized_fe_functions

    def _localize_non_fe_functions(self, mesh_entity):
        for k in self.non_fe_functions.keys():
            self._localized_non_fe_functions[k] = self.non_fe_functions[k].localize(mesh_entity)
        return self._localized_non_fe_functions

    def set_mesh(self, mesh, subdomain=None):
        self._mesh = mesh
        self._subdomain = subdomain
        for fef in self.fe_functions.values():
            if hasattr(fef, "set_mesh"):
                self.mapping = fef.set_mesh(mesh, subdomain)

    @abc.abstractmethod
    def implements_quadrature_on(self, entity_type=None):
        """
        Implementation should give information whether a cell_linear_form etc is provided or not.
        :param entity_type: Form.cells, .interior_faces, .exterior_faces to query if such integrals are computed.
        In this case 'cells' refers to the top level entities (in 1d these are lines); 'faces' denotes
        the second highes entities.
        """
        raise Exception("Abstract method called!")

    def mesh_entity_iterator(self, topological_dim=None):
        return self._mesh.get_mesh_entities(topological_dim=topological_dim)


class Functional(Form):
    """
    Represesent functionals, usually in the sense of an integral of some functions over a given domain.
    This class is an abstract class! However, only the local (element / inter-element) computations have to be
    implemented. Have a look at the types
    `CellEvalDataBase`, `InteriorFaceEvalDataBase` and `ExteriorFaceEvalDataBase`
    for information that is available for these local contributions.
    In case you are sure to do e.g, only cell-evaluations but no cell-boundary evaluations you may
    omit the implementation of the "unneeded" methods via a "pass" or raising an NotImplementedError().
    """
    def __init__(self, quadrature, mapping=None, fe_functions=None, non_fe_functions=None, mesh=None, subdomain=None):
        Form.__init__(self, quadrature, fe_functions=fe_functions, non_fe_functions=non_fe_functions, mapping=mapping,
                      mesh=mesh, subdomain=subdomain)
        if self._mesh is None:
            for fef in self.fe_functions.values():
                if hasattr(fef, "get_mesh") and hasattr(fef, "get_subdomain"):
                    try:
                        self.set_mesh(fef.get_mesh(), fef.get_subdomain())
                        break
                    except TypeError:
                        pass

    @abc.abstractmethod
    def local_cell_functional(self, cell_eval_data_functional):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_interior_face_functional(self, interior_face_eval_data_functional):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_exterior_face_functional(self, exterior_face_eval_data_functional):
        raise Exception("Abstract method called!")

    def get_cell_eval_data_functional(self, mesh_entity, params=None, mapping=None):
        local_mapping = None
        if mapping is None:
            mapping = self.mapping
        if mapping is not None:
            local_mapping = mapping.localize(mesh_entity)

        return CellEvalDataBase(self._localize_fe_functions(mesh_entity),
                                self._localize_non_fe_functions(mesh_entity),
                                local_mapping,
                                self.quadrature,
                                params)

    def get_interior_face_eval_data_functional(self, mesh_entities, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        local_mappings = None
        if mapping is None:
            mapping = self.mapping
        if mapping is not None:
            local_mappings = [mapping.localize(mesh_entities[0]),
                              mapping.localize(mesh_entities[1])]

        return InteriorFaceEvalDataBase(
            [self._localize_fe_functions(mesh_entities[0]),
             self._localize_fe_functions(mesh_entities[1])],
            [self._localize_non_fe_functions(mesh_entities[0]),
             self._localize_non_fe_functions(mesh_entities[1])],
            local_mappings,
            self.quadrature,
            params)

    def get_exterior_face_eval_data_functional(self, mesh_entity, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        local_mapping = None
        if mapping is None:
            mapping = self.mapping
        if mapping is not None:
            local_mapping = mapping.localize(mesh_entity)
        return ExteriorFaceEvalDataBase(self._localize_fe_functions(mesh_entity),
                                        self._localize_non_fe_functions(mesh_entity),
                                        local_mapping,
                                        self.quadrature,
                                        params)


class LinearForm(Form):
    """
    Represesents linear forms (in the FEM sense),over a given domain.
    This class is an abstract class! However, only the local (element / inter-element) computations have to be
    implemented. Have a look at the types
    `CellEvalDataLinearForm`, `InteriorFaceEvalDataLinearForm` and `ExteriorFaceEvalDataLinearForm`
    for information that is available for these local contributions.
    """
    def __init__(self, test_function_space, quadrature, fe_functions=None, non_fe_functions=None, mapping=None,
                 mesh=None, subdomain=None):
        Form.__init__(self, quadrature, mapping=mapping, fe_functions=fe_functions, non_fe_functions=non_fe_functions)
        self.test_function_space = test_function_space
        if self.mapping is None:
            self.mapping = self.test_function_space.get_mapping()
        if mesh is not None:
            self.test_function_space.set_mesh(mesh, subdomain=subdomain)
        else:
            self.set_mesh(self.test_function_space.get_mesh(), self.test_function_space.get_subdomain())

    @abc.abstractmethod
    def local_cell_linear_form(self, cell_eval_data_linear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_interior_face_linear_form(self, interior_face_eval_data_linear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_exterior_face_linear_form(self, exterior_face_eval_data_linear_form):
        raise Exception("Abstract method called!")

    def set_mesh(self, mesh, subdomain=None):
        Form.set_mesh(self, mesh, subdomain)
        self.test_function_space.set_mesh(mesh, subdomain)

    def get_cell_eval_data_linear_form(self, mesh_entity, params=None, mapping=None):
        local_mapping = None
        if mapping is None:
            mapping = self.mapping
        if mapping is not None:
            local_mapping = mapping.localize(mesh_entity)
        return CellEvalDataLinearForm(self.test_function_space.localize(mesh_entity),
                                      self._localize_fe_functions(mesh_entity),
                                      self._localize_non_fe_functions(mesh_entity),
                                      local_mapping,
                                      self.quadrature,
                                      params)

    def get_interior_face_eval_data_linear_form(self, mesh_entities, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        local_mappings = None
        if mapping is None:
            mapping = self.mapping
            if mapping is not None:
                local_mappings = [mapping.localize(mesh_entities[0]),
                                  mapping.localize(mesh_entities[1])]
        return InteriorFaceEvalDataLinearForm(
            [self.test_function_space.localize(mesh_entities[0]),
             self.test_function_space.localize(mesh_entities[1])],
            [self._localize_fe_functions(mesh_entities[0]),
             self._localize_fe_functions(mesh_entities[1])],
            [self._localize_non_fe_functions(mesh_entities[0]),
             self._localize_non_fe_functions(mesh_entities[1])],
            local_mappings,
            self.quadrature,
            params)

    def get_exterior_face_eval_data_linear_form(self, mesh_entity, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        local_mapping = None
        if mapping is None:
            mapping = self.mapping
            if mapping is not None:
                local_mapping = mapping.localize(mesh_entity)
        return ExteriorFaceEvalDataLinearForm(self.test_function_space.localize(mesh_entity),
                                              self._localize_fe_functions(mesh_entity),
                                              self._localize_non_fe_functions(mesh_entity),
                                              local_mapping,
                                              self.quadrature,
                                              params)

    def get_local_size(self, mesh_entity):
        return self.test_function_space.get_number_of_global_element_dofs(mesh_entity)

    def get_test_function_space_dim(self):
        return self.test_function_space.function_space_dim()


class BilinearForm(Form):
    """
    Represesents bilinear forms (in the FEM sense),over a given domain.
    This class is an abstract class! However, only the local (element / inter-element) computations have to be
    implemented. Have a look at the types
    `CellEvalDataBilinearForm`, `InteriorFaceEvalDataBilinearForm` and `ExteriorFaceEvalDataBilinearForm`
    for information that is available for these local contributions.
    """
    def __init__(self, test_function_space, trial_function_space, quadrature, mapping=None, fe_functions=None,
                 non_fe_functions=None, mesh=None, subdomain=None):
        Form.__init__(self, quadrature, mapping=mapping, fe_functions=fe_functions, non_fe_functions=non_fe_functions)
        self.test_function_space = test_function_space
        self.trial_function_space = trial_function_space
        if self.mapping is None:
            self.mapping = self.test_function_space.get_mapping()
        if mesh is not None:
            self.test_function_space.set_mesh(mesh, subdomain=subdomain)
            self.trial_function_space.set_mesh(mesh, subdomain=subdomain)
        else:
            self.set_mesh(self.test_function_space.get_mesh(), self.test_function_space.get_subdomain())

    @abc.abstractmethod
    def local_cell_bilinear_form(self, cell_eval_data_bilinear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_interior_face_bilinear_form(self, interior_face_eval_data_bilinear_form):
        raise Exception("Abstract method called!")

    @abc.abstractmethod
    def local_exterior_face_bilinear_form(self, exterior_face_eval_data_bilinear_form):
        raise Exception("Abstract method called!")

    def set_mesh(self, mesh, subdomain=None):
        Form.set_mesh(self, mesh, subdomain)
        self.test_function_space.set_mesh(mesh, subdomain)
        self.trial_function_space.set_mesh(mesh, subdomain)

    def get_cell_eval_data_bilinear_form(self, mesh_entity, params=None, mapping=None):
        local_mapping = None
        if mapping is None:
            mapping = self.mapping
            if mapping is not None:
                local_mapping = mapping.localize(mesh_entity)
        return CellEvalDataBilinearForm(self.test_function_space.localize(mesh_entity),
                                        self.trial_function_space.localize(mesh_entity),
                                        self._localize_fe_functions(mesh_entity),
                                        self._localize_non_fe_functions(mesh_entity),
                                        local_mapping,
                                        self.quadrature,
                                        params)

    def get_interior_face_eval_data_bilinear_form(self, mesh_entities, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        local_mappings = None
        if mapping is None:
            mapping = self.mapping
            if mapping is not None:
                local_mappings = [mapping.localize(mesh_entities[0]),
                                  mapping.localize(mesh_entities[1])]
        return InteriorFaceEvalDataBilinearForm(
            [self.test_function_space.localize(mesh_entities[0]),
             self.test_function_space.localize(mesh_entities[1])],
            [self.trial_function_space.localize(mesh_entities[0]),
             self.trial_function_space.localize(mesh_entities[1])],
            [self._localize_fe_functions(mesh_entities[0]),
             self._localize_fe_functions(mesh_entities[1])],
            [self._localize_non_fe_functions(mesh_entities[0]),
             self._localize_non_fe_functions(mesh_entities[1])],
            local_mappings,
            self.quadrature,
            params)

    def get_exterior_face_eval_data_bilinear_form(self, mesh_entity, params=None, mapping=None):
        # TODO: reason about "selecting" the exterior face
        # should there be a list of boundary indicators to be respected?
        # should the boundary indicators be respected by the implementations of the (bi)linear forms?
        # should the mesh entity be the actual face and not the corresponding cell?
        # should the number of the faces be given?
        # should this be handled in the cell_*-stuff?
        # ...
        local_mapping = None
        if mapping is None:
            mapping = self.mapping
            if mapping is not None:
                local_mapping = mapping.localize(mesh_entity)
        return ExteriorFaceEvalDataBilinearForm(self.test_function_space.localize(mesh_entity),
                                                self.trial_function_space.localize(mesh_entity),
                                                self._localize_fe_functions(mesh_entity),
                                                self._localize_non_fe_functions(mesh_entity),
                                                local_mapping,
                                                self.quadrature,
                                                params)

    def get_local_shape(self, mesh_entity):
        return (self.test_function_space.get_number_of_global_element_dofs(mesh_entity),
                self.trial_function_space.get_number_of_global_element_dofs(mesh_entity))

    def get_test_function_space_dim(self):
        return self.test_function_space.function_space_dim()

    def get_trial_function_space_dim(self):
        return self.trial_function_space.function_space_dim()


class FormCollection(object):
    """
    This class represents a collection of forms by means of iterators. It is possible to add forms to and remove
    forms from the collections. Access to elements is only provided via iterators over dict-values.
    These iterators are fitlered iterators and are implemented for functionals, linear forms and bilinear forms.

    The only purpose of this class is to aggregate forms of different kind in a siingle structure and provide
    iterators of them to e.g. an assembly algortihm.
    """
    def __init__(self, form_dict=None):
        self._form_dict = form_dict
        if form_dict is None:
            self._form_dict = {}

    def set_form(self, key, form):
        """
        Register the given form under key.
        :param key: almost anything, but it is a good idea to provide an integer as index or a "name" of the form
        :param form: something that inherits from Form.
        :return: None
        """
        self._form_dict[key] = form

    def unset_form(self, key):
        """
        Removes (and returns) the form stored under key from the collection.
        :param key: the key of the object that is to be removed.
        :return: the object stored under key
        """
        self._form_dict.pop(key)

    def functional_iterator(self):
        return filter(lambda f: hasattr(f, "local_cell_functional"), self._form_dict.values())

    def linear_form_iterator(self):
        return filter(lambda f: hasattr(f, "local_cell_linear_form"), self._form_dict.values())

    def bilinear_form_iterator(self):
        return filter(lambda f: hasattr(f, "local_cell_bilinear_form"), self._form_dict.values())
