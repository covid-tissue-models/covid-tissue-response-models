"""
This module defines the CompuCell3D base steppable class for the framework.

Written by T.J. Sego, Ph.D.
"""

from collections import namedtuple
from typing import *

from cc3d.core.PySteppables import *
from cc3d.cpp import CompuCell

ode_manager_key = 'std_ode_manager_steppable'
_registry_key = 'nCoVSteppableBase_registry'


class NoODEManagerException(Exception):
    pass


class nCoVSteppableBase(SteppableBasePy):
    """
    The base class of all steppables in the viral infection framework.

    This is an intermediate base class between SteppableBasePy, which is the base class used in model implementation
    in the typical cc3d fashion, and a steppable-based model implementation in this framework.
    This class defines methods and features that support basic model implementation via a steppable without complete
    knowledge or control of how, or in what context, the steppable will be used.
    For example, in a typical cc3d model implementation, a modeler can define and call routines associated with the
    creation of a new cell (e.g., attach an attribute in the cell dictionary). However, a steppable derived from this
    class defines routines associated with the creation of a new cell independently of what steppable creates the cell.

    Every derived class must define the class attribute *unique_key*, which makes it accessible by all other derived
    steppable instances when registered in the shared steppable dictionary.
    All registered steppables derived from this class are available in *shared_steppable_vars* to all steppables at,
    and after, cc3d calls the steppable method *start*.

    Defines attributes ``step_period`` and ``voxel_length`` for unit conversion of time and length.
    These attributes are connected to class variables, so their values are the same across all instances of
    derived classes.
    Derived classes should convert their parameters to cc3d units (`e.g.`, from minutes to steps) using these.

    Defines a callback *on_new_cell*, which is called on every registered derived class when any derived class
    creates a new cell in the typical cc3d fashion using *new_cell*.
    Derived classes override this method to implement routines associated with the creation of a new cell.
    Calls are made in order of registration with cc3d.

    Defines a callback *on_set_cell*, which is called on every registered derived class when any derived class
    sets the type of a cell to a new value using *set_cell_type* (rather than the typical cc3d fasion of setting the
    attribute *type* on a cell instance).
    As such, all derived class should use *set_cell_type* to broadcast changes in cell type to all registered
    derived classes of a simulation.
    Derived classes override this method to implement routines associated with the change in a cell type.
    Calls are made in order of registration with cc3d.

    Defines an interface to an ODE manager, which automatically instantiates and destroys instances of ODE models
    attached to cells when they are created or their type is changed.
    ODE models, whether free-floating or attached to cells of a specific type, should be registered using
    *register_ode_model*.
    A steppable that defines and registers an ODE model, whether free-floating or attached to cells, should perform
    stepping of the ODE model using *timestep_ode_model*, rather than the typical cc3d fashion of one steppable
    performing stepping all ODE models using *timestep_sbml*.
    """

    unique_key = None
    """
    Unique key to find a steppable in the shared steppable dictionary.
    
    This must be set by a derived class so that other steppables can find it if registered with CC3D
    """

    _step_period = 60.0

    @property
    def step_period(self) -> float:
        """
        Unit conversion for time, in units seconds per step
        """
        return nCoVSteppableBase._step_period

    @step_period.setter
    def step_period(self, _val: float):
        if _val <= 0:
            raise ValueError('Step period must be positive')
        nCoVSteppableBase._step_period = _val

    _voxel_length = 1.0

    @property
    def voxel_length(self) -> float:
        """
        Unit conversion for space, in micrometers per voxel side length
        """
        return nCoVSteppableBase._voxel_length

    @voxel_length.setter
    def voxel_length(self, _val: float):
        if _val <= 0:
            raise ValueError('Site length must be positive')
        nCoVSteppableBase._voxel_length = _val

    def core_init(self, reinitialize_cell_types=True):
        super().core_init(reinitialize_cell_types)

        # Share self
        if self.unique_key is None:
            raise AttributeError("Steppable key not set")
        elif self.unique_key in self.shared_steppable_vars.keys():
            # Handling the occurence where core_init is called multiple times on consecutive runs for some steppables
            return
        else:
            print("Registering nCoVSteppableBase:", self.unique_key)
        self.shared_steppable_vars[self.unique_key] = self
        self._register_base()

        # Register ODEManagerSteppable if necessary
        from cc3d import CompuCellSetup
        steppable_registry = CompuCellSetup.persistent_globals.steppable_registry
        if ODEManagerSteppable.__name__ not in steppable_registry.steppableDict.keys():
            CompuCellSetup.register_steppable(ODEManagerSteppable(frequency=1))

    @property
    def _base_registry(self):
        if _registry_key not in self.shared_steppable_vars.keys():
            self.shared_steppable_vars[_registry_key] = []
        return self.shared_steppable_vars[_registry_key]

    def _register_base(self):
        if self not in self._base_registry:
            self._base_registry.append(self)

    def _unregister_base(self):
        if self in self._base_registry:
            self._base_registry.remove(self)

    def new_cell(self, cell_type=0) -> CompuCell.CellG:
        new_cell = super().new_cell(cell_type)

        if self.ode_manager is not None:
            for model_entry in self.ode_manager.models_by_cell_type(self.get_type_name_by_cell(new_cell)):
                self._new_cell_ode_model(cell=new_cell, model_name=model_entry.model_name)

        # Do callback on all registered subclasses
        # If a subclass instance returns False, try again in a subsequent iteration
        steppables = self._base_registry.copy()
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

    def set_cell_type(self, cell: CompuCell.CellG, _type_id: int) -> None:
        """
        Sets the type of a cell.

        This is a supplement to the typical cc3d fasion of setting the type of a cell using *cell.type*.

        A subsequent call to *on_set_cell_type* is made on all registered derived classes.

        :param cell: a cell
        :param _type_id: the new type id of the cell
        :return: None
        """
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
        steppables = self._base_registry.copy()
        while steppables:
            steppables = [x for x in steppables if x.on_set_cell_type(cell, old_type_id) is False]

    def on_set_cell_type(self, cell: CompuCell.CellG, old_type: int) -> Union[None, bool]:
        """
        A callback on every occasion that a subclass calls *set_cell_type*.

        If the returned value is False, then the callback is deferred and will be called again in a subsequent
        iteration over all subclasses registered with cc3d.
        The callback is called after the cell type has been changed and all ODE models have been updated according to
        the ode manager.

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
        """
        Reference to the ODEManagerSteppable instance of a simulation.
        """
        if ode_manager_key in self.shared_steppable_vars.keys():
            return self.shared_steppable_vars[ode_manager_key]
        return None

    @property
    def ode_model_names(self) -> List[str]:
        """
        List of registered ode models
        """
        if self.ode_manager is not None:
            return self.ode_manager.mode_names
        raise NoODEManagerException

    def ode_models_by_cell_type(self, _cell_type: str):
        """
        List of registered ode model names associated with a cell type name in a simulation.
        """
        if self.ode_manager is not None:
            return self.ode_manager.models_by_cell_type(_cell_type)
        raise NoODEManagerException

    def cell_types_by_ode_model(self, _model_name: str) -> Union[None, List[int]]:
        """
        List of cell type ids associated with a registered ode model in a simulation.
        """
        if self.ode_manager is not None:
            return self.ode_manager.cell_types_by_model(_model_name)
        raise NoODEManagerException

    def register_ode_model(self,
                           model_name: str,
                           model_fcn: Callable,
                           ics_fcn: Callable = None,
                           cell_types: Union[str, Iterable[str], None] = None,
                           model_type: str = 'antimony',
                           step_size: float = 1.0) -> None:
        """
        Registers an ode model with the framework and cc3d.

        This can be used to register a free-floating ode model, or an ode model attached to cells of a set of types.

        Free-floating ode models are automatically instantiated on registration.

        Ode models attached to cells are automatically managed by the framework by their type.
        *E.g.*, if an ode model is registered for cell type "B" and a cell of type "A" is changed to type "B" using
        *set_cell_type*, then the framework automatically instantiates an instance of the ode model and attaches it to
        the cell.

        :param model_name: name of the ode model
        :param model_fcn: model string generator function.
            Free-floating ode model generator functions take no arguments.
            Generator functions for ode models attached to cells take the argument of the cell to which the ode model is
            being attached.
            Generator functions return a multiline string of the ode model in antimony, cellml or sbml model syntax.
        :param ics_fcn: initial conditions generator function (optional).
            Free-floating ode model initial conditions generator functions take no arguments.
            Generator functions for ode model initial conditions attached to cells take the argument of the cell to
            which the ode model is being attached.
            Generator functions return a dictionary with key: value pairs of sbml model variable string name: variable
            value.
        :param cell_types: corresponding name of cell type or list of names of cell types if attached to cells;
            if not specified, the ode model is registered as free-flaoting.
        :param model_type: string name of model language ('antimony', 'cellml' or 'sbml'); default 'antimony'
        :param step_size: step size of ode model; default 1.0.
        :return: None
        """
        if self.ode_manager is not None:
            return self.ode_manager.register_model(model_name=model_name,
                                                   model_fcn=model_fcn,
                                                   ics_fcn=ics_fcn,
                                                   cell_types=cell_types,
                                                   model_type=model_type,
                                                   step_size=step_size)
        raise NoODEManagerException

    def timestep_ode_model(self, model_name: str) -> None:
        """
        Integrate an ode model one step in time.

        :param model_name: name of the ode model
        :return: None
        """
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
        """
        List of registered ode models
        """
        return [x for x in self._ode_models.keys()]

    def models_by_cell_type(self, _cell_type: str) -> List[ODEModelEntry]:
        """
        List of registered ode model names associated with a cell type in a simulation
        """
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
        """
        List of cell type ids associated with a registered ode model in a simulation.
        """
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
        """
        Registers an ode model with the framework and cc3d.

        This can be used to register a free-floating ode model, or an ode model attached to cells of a set of types.

        Free-floating ode models are automatically instantiated on registration.

        Ode models attached to cells are automatically managed by the framework by their type.
        *E.g.*, if an ode model is registered for cell type "B" and a cell of type "A" is changed to type "B" using
        *set_cell_type*, then the framework automatically instantiates an instance of the ode model and attaches it to
        the cell.

        :param model_name: name of the ode model
        :param model_fcn: model string generator function.
            Free-floating ode model generator functions take no arguments.
            Generator functions for ode models attached to cells take the argument of the cell to which the ode model is
            being attached.
            Generator functions return a multiline string of the ode model in antimony, cellml or sbml model syntax.
        :param ics_fcn: initial conditions generator function (optional).
            Free-floating ode model initial conditions generator functions take no arguments.
            Generator functions for ode model initial conditions attached to cells take the argument of the cell to
            which the ode model is being attached.
            Generator functions return a dictionary with key: value pairs of sbml model variable string name: variable
            value.
        :param cell_types: corresponding name of cell type or list of names of cell types if attached to cells;
            if not specified, the ode model is registered as free-flaoting.
        :param model_type: string name of model language ('antimony', 'cellml' or 'sbml'); default 'antimony'
        :param step_size: step size of ode model; default 1.0.
        :return: None
        """
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
        """
        Integrate an ode model one step in time.

        :param model_name: name of the ode model
        :return: None
        """
        if model_name not in self.model_names:
            raise ValueError('Model not registered:', model_name)
        o = self._ode_models[model_name]
        if o.cell_types is None:
            from cc3d.CompuCellSetup import persistent_globals as pg
            pg.free_floating_sbml_simulators[model_name].timestep()
        else:
            for cell in self.cell_list_by_type(*self.cell_types_by_model(model_name)):
                dict_attrib = CompuCell.getPyAttrib(cell)
                dict_attrib['SBMLSolver'][model_name].timestep()
