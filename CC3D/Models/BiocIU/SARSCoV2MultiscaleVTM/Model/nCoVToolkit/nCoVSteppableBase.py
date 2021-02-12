# This is a general library of CompuCell3D steppable classes for the shared coronavirus modeling and simulation project
# hosted by the Biocomplexity Institute at Indiana University

from collections import namedtuple
from typing import *

from cc3d.core.PySteppables import *
from cc3d.cpp import CompuCell

ode_manager_key = 'std_ode_manager_steppable'


class NoODEManagerException(Exception):
    pass


class nCoVSteppableBase(SteppableBasePy):

    unique_key = None
    """
    Unique key to find a steppable in the shared steppable dictionary.
    
    This must be set by a derived class so that other steppables can find it if registered with CC3D
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def core_init(self, reinitialize_cell_types=True):
        super().core_init(reinitialize_cell_types)

        # Share self
        if self.unique_key is None:
            raise AttributeError("Steppable key not set")
        elif self.unique_key in self.shared_steppable_vars.keys():
            # Currently we can't enforce this, since shared_steppable_vars does not reset at every simulation instance
            # raise ValueError(f"Steppable key is not unique: {self.unique_key}")
            pass
        else:
            print("Registering nCoVSteppableBase:", self.unique_key)
        self.shared_steppable_vars[self.unique_key] = self

        # Register ODEManagerSteppable if necessary
        from cc3d.CompuCellSetup import persistent_globals
        steppable_registry = persistent_globals.steppable_registry
        if not any([isinstance(x, ODEManagerSteppable) for x in steppable_registry.allSteppables()]):
            steppable_registry.registerSteppable(ODEManagerSteppable(frequency=1))

    def step(self, mcs):
        pass

    def finish(self):
        pass

    def new_cell(self, cell_type=0):
        new_cell = super().new_cell(cell_type)

        if self.ode_manager is not None:
            for model_entry in self.ode_manager.models_by_cell_type(self.get_type_name_by_cell(new_cell)):
                self._new_cell_ode_model(cell=new_cell, model_name=model_entry.model_name)

        # Do callback on all registered subclasses
        # If a subclass instance returns False, try again in a subsequent iteration
        from cc3d.CompuCellSetup import persistent_globals
        steppables = [x for x in persistent_globals.steppable_registry.allSteppables()
                      if isinstance(x, nCoVSteppableBase)]
        while steppables:
            steppables = [x for x in steppables if x.on_new_cell(new_cell) is False]

        return new_cell

    def on_new_cell(self, _new_cell: CompuCell.CellG) -> Union[None, bool]:
        """
        A callback on every occasion that a subclass calls new_cell. If the returned value is False, then the callback
        is deferred and will be called again in a subsequent iteration over all subclasses registered with cc3d.

        :param _new_cell: newly created cell
        :return: None if not deferring; False if deferring
        """
        pass

    def set_cell_type(self, cell: CompuCell.CellG, _type_id: int):
        if cell.type == _type_id:
            return

        # Remove ODE models that are no longer relevant to this cell
        if self.ode_manager is not None:
            for model_entry in self.ode_manager.models_by_cell_type(self.get_type_name_by_cell(cell)):
                cell_type_ids = [getattr(self, ct.upper()) for ct in model_entry.cell_types]
                if _type_id not in cell_type_ids:
                    self.delete_sbml_from_cell(model_name=model_entry.model_name,
                                               cell=cell)

        old_type_id = cell.type
        cell.type = _type_id

        # Add ODE models that are now relevant to this cell
        if self.ode_manager is not None:
            for model_entry in self.ode_manager.models_by_cell_type(self.get_type_name_by_cell(cell)):
                if not hasattr(cell.sbml, model_entry.model_name):
                    self._new_cell_ode_model(cell=cell, model_name=model_entry.model_name)

        # Do callback on all registered subclasses
        # If a subclass instance returns False, try again in a subsequent iteration
        from cc3d.CompuCellSetup import persistent_globals
        steppables = [x for x in persistent_globals.steppable_registry.allSteppables()
                      if isinstance(x, nCoVSteppableBase)]
        while steppables:
            steppables = [x for x in steppables if x.on_set_cell_type(cell, old_type_id) is False]

    def on_set_cell_type(self, cell: CompuCell.CellG, old_type: int) -> Union[None, bool]:
        """
        A callback on every occasion that a subclass calls set_cell_type. If the returned value is False, then the
        callback is deferred and will be called again in a subsequent iteration over all subclasses registered with
        cc3d. The callback is called after the cell type has been changed and all ODE models have been updated
        according to the ode manager

        :param cell: a cell with a newly assigned type
        :param old_type: previous type of the cell
        :return: None if not deferring; False if deferring
        """
        pass

    def _new_cell_ode_model(self, cell: CompuCell.CellG, model_name: str):
        if self.ode_manager is None:
            raise NoODEManagerException
        model_entry = self.ode_manager._ode_models[model_name]
        initial_conditions = None
        if model_entry.ics_fcn is not None:
            initial_conditions = model_entry.ics_fcn(cell)
        if model_entry.model_type == 'antimony':
            loader = self.add_antimony_to_cell
        elif model_entry.model_name == 'sbml':
            loader = self.add_sbml_to_cell
        elif model_entry.model_name == 'cellml':
            loader = self.add_cellml_to_cell
        else:
            raise TypeError('Valid model types are antimony, sbml and cellml')
        loader(model_string=model_entry.model_fcn(cell),
               model_name=model_entry.model_name,
               cell=cell,
               step_size=model_entry.step_size,
               initial_conditions=initial_conditions)

    @property
    def ode_manager(self):
        if ode_manager_key in self.shared_steppable_vars.keys():
            return self.shared_steppable_vars[ode_manager_key]
        return None

    def ode_models_by_cell_type(self, _cell_type: str):
        if self.ode_manager is not None:
            return self.ode_manager.models_by_cell_type(_cell_type)
        raise NoODEManagerException

    def cell_types_by_ode_model(self, _model_name: str) -> Union[None, List[int]]:
        if self.ode_manager is not None:
            return self.ode_manager.cell_types_by_model(_model_name)
        raise NoODEManagerException

    def register_ode_model(self,
                           model_name: str,
                           model_fcn: Callable,
                           ics_fcn: Callable = None,
                           cell_types: Union[str, Iterable[str], None] = None,
                           model_type: str = 'antimony',
                           step_size: float = 1.0):
        if self.ode_manager is not None:
            return self.ode_manager.register_model(model_name=model_name,
                                                   model_fcn=model_fcn,
                                                   ics_fcn=ics_fcn,
                                                   cell_types=cell_types,
                                                   model_type=model_type,
                                                   step_size=step_size)
        raise NoODEManagerException

    def timestep_ode_model(self, model_name: str):
        if self.ode_manager is None:
            raise NoODEManagerException
        self.ode_manager.timestep_model(model_name=model_name)


ODEModelEntry = namedtuple(typename='ODEModelEntry',
                           field_names=['model_name', 'model_fcn', 'ics_fcn', 'cell_types', 'model_type', 'step_size'])


class ODEManagerSteppable(nCoVSteppableBase):
    """
    Hub for accessible free-floating and cell ODE models

    This steppable does not step an ODE model; stepping an ODE model is the responsibility of its owner

    This steppable contains a registry for registering free-floating and cell ODE models, which is used by
    nCoVSteppableBase to automatically manage ODE models (*e.g.*, when creating a cell or changing its type) and make
    working with them more robust to support an arbitrary array of steppables

    This steppable is automatically registered, and does not need to be registered manually
    """

    unique_key = ode_manager_key

    def __init__(self, frequency=1):
        nCoVSteppableBase.__init__(self, frequency)

        self._ode_models: Dict[str, ODEModelEntry] = {}

    @property
    def model_names(self) -> List[str]:
        return [x for x in self._ode_models.keys()]

    def models_by_cell_type(self, _cell_type: str) -> List[ODEModelEntry]:
        x = []
        for y in self._ode_models.values():
            cell_types = y.cell_types
            if cell_types is not None:
                if isinstance(cell_types, str) and cell_types == _cell_type:
                    x.append(y)
                elif isinstance(cell_types, Iterable) and _cell_type in cell_types:
                    x.append(y)
        return x

    def cell_types_by_model(self, _model_name: str) -> Union[None, List[int]]:
        ode_model_entry = self._ode_models[_model_name]
        if ode_model_entry.cell_types is None:
            return None
        if isinstance(ode_model_entry.cell_types, str):
            return [getattr(self, ode_model_entry.cell_types.upper())]
        return [getattr(self, x.upper()) for x in self._ode_models[_model_name].cell_types]

    def register_model(self,
                       model_name: str,
                       model_fcn: Callable,
                       ics_fcn: Callable = None,
                       cell_types: Union[str, Iterable[str], None] = None,
                       model_type: str = 'antimony',
                       step_size: float = 1.0):
        if model_type.lower() not in ['antimony', 'sbml', 'cellml']:
            raise TypeError('Valid model types are antimony, sbml and cellml')
        elif model_name in self.model_names:
            raise AttributeError('Model name is not unique')
        elif cell_types is not None:
            err_str = 'Cell type name not registered:'
            if isinstance(cell_types, str) and not hasattr(self, cell_types.upper()):
                raise ValueError(err_str, cell_types)
            elif isinstance(cell_types, Iterable):
                for ct in cell_types:
                    if not hasattr(self, ct.upper()):
                        raise ValueError(err_str, ct)

        ode_model_entry = ODEModelEntry(model_name=model_name,
                                        model_fcn=model_fcn,
                                        ics_fcn=ics_fcn,
                                        cell_types=cell_types,
                                        model_type=model_type.lower(),
                                        step_size=step_size)

        # Automatically add free-floating models
        if cell_types is None:
            initial_conditions = None
            if ics_fcn is not None:
                initial_conditions = ics_fcn()
            loader = {'antimony': self.add_free_floating_antimony,
                      'sbml': self.add_free_floating_sbml,
                      'cellml': self.add_free_floating_cellml}[model_type]
            loader(model_string=model_fcn(),
                   model_name=model_name,
                   step_size=step_size,
                   initial_conditions=initial_conditions)

        self._ode_models[model_name] = ode_model_entry

    def timestep_model(self, model_name: str):
        if model_name not in self.model_names:
            raise ValueError('Model not registered:', model_name)
        o = self._ode_models[model_name]
        if o.cell_types is None:
            from cc3d.CompuCellSetup import persistent_globals as pg
            pg.free_floating_sbml_simulators[model_name].timestep()
        else:
            for cell in self.cell_list_by_type(*o.cell_types):
                dict_attrib = CompuCell.getPyAttrib(cell)
                dict_attrib['SBMLSolver'].timestep()
