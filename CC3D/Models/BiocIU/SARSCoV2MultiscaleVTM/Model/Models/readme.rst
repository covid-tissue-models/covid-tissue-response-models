Usage
=====

This framework functions as a typical CompuCell3D (CC3D) simulation, but with support for a library of shared,
community-developed model modules.
In a typical CC3D model implementation, a user defines python steppables that describe various aspects of a model and
registers them with CC3D, in which case they are simulated in tandem with built-in plugins and steppables as also
specified by the user. In addition to typical CC3D simulation capabilities, this framework also provides a library
of shared, community-developed model modules that the user can selectively utilize, manipulate, and modify according
to their own modeling project.

Using a steppable defined in a module in this collection is as simple as using any other steppable in CC3D.
To use a steppable in simulation defined in a module in this collection, simply import from the module
``Models`` and add it to the loading procedures in ``Simulation/ViralInfectionVTM.py``. For example,
to use the steppable ``SimpleRecoverySteppable`` defined in the model module ``Models.RecoverySimple``,
add the following to ``Simulation/ViralInfectionVTM.py``,

.. code-block:: python

    from Models.RecoverySimple.RecoverySteppables import SimpleRecoverySteppable
    CompuCellSetup.register_steppable(steppable=SimpleRecoverySteppable(frequency=1))

Comparing this procudure to those included in ``Simulation/ViralInfectionVTM.py``, we see that not much is
different. The only difference is specifying the location of the Python script containing the steppable that
we would like to use. This way, two model modules can define steppables in Python scripts of the same name
without overwriting each other (`e.g.`, ``Models.RecoverySimple.RecoverySteppables`` vs.
``Models.RecoveryNeighbor.RecoverySteppables``). The only necessarily unique aspect of a particular model
module is the name of its containing directory (*e.g.*, ``Models.RecoverySimple`` vs.
``Models.RecoveryNeighbor``). This scheme isolates model-specific development to the directory in which
the add-on model is defined, and modularizes the overall simulation framework into `shareable`,
`interchangable` model components.

For available features, options and other information about a particular model module, please refer to the
documentation provided within the module.

In a typical CC3D model implementation, a user defines a python steppable by deriving from the CC3D class
`SteppableBasePy`, or one of its derived classes (`e.g.`, `SecretionBasePy`).
For example, a typical specification of a steppable called `MySteppable` would look like the following,

.. code-block:: python

    from cc3d.core.PySteppables import *

    class MySteppable(SteppableBasePy):
        def __init__(self, frequency=1):
            SteppableBasePy.__init__(self, frequency)

To exploit the full functionality of this framework, substitute the usage of ``SteppableBasePy`` with
``nCoVSteppableBase`` (named after the pandemic that catalyzed this project), which is an intermediate class on top of
``SteppableBasePy`` that has all of the functionality of ``SteppableBasePy`` plus built-in functionality for this
framework. This can be accomplished as follows, using the previous example of defining ``MySteppable``,

.. code-block:: python

    from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase

    class MySteppable(nCoVSteppableBase):
        def __init__(self, frequency=1):
            nCoVSteppableBase.__init__(self, frequency)

This functionality will work out-of-the-box if you use the provided script ``ViralInfectionVTM.py`` in the directory
``Simulation`` for importing and registering steppables with CC3D. Specifically, ``ViralInfectionVTM.py`` issues a few
requisite commands that make the framework accessible before importing ``CompuCellSetup`` from ``cc3d``. So if you
would like to use a completely different script, or write one `de novo`, be sure to issue the same commands.

``nCoVSteppableBase`` is a special class that allows you to exploit both the basic functionality of CC3D `and` the
functionality of the shared models that are distributed with this framework.
While each module steppable can define ways in which you can interact with it, modify its internal data, or connect it
with other modules, some information about a simulation is intrinsically shared among many, if not all, steppables of
a simulation.
For example, all CC3D simulations work with the same inventory of cells, each of which can have only one type among
a set of cell types in the simulation according to your simulation specification.
This commonality motivates why every class that uses ``SteppableBasePy`` refers to the same cell inventory when it
accesses its attribute ``cell_list`` (`i.e.`, two steppables may define different operations on the cell inventory, but
each is still performing their operations on the `same` cell inventory).
``nCoVSteppableBase`` builds upon this same principle in important ways so that you can issue framework-wide changes
and events that are reflected in whatever steppables you've incorporated into your simulation, whether they are
steppables you have designed, or they are distributed in the shared library. Issuing and handling events are both
described later in this text.

CC3D uses space and time units of `voxel length` and `step`, respectively, which makes a CC3D simulation more
robust in that the size of a voxel or the period of a simulation step are whatever you decide (if you care at all).
For example, if I decide that the voxels of my simulation are 1 x 1 x 1 microns and a step is one minute, then
a lattice of 1000 x 1000 x 1000 voxels and 60 steps represent one hour in a 1 x 1 x 1 mm volume.
For any model specification with model parameters, a steppable must then convert their model parameters to those of
CC3D.
``nCoVSteppableBase`` provides shared data with all steppables for making their unit conversions, and so that
you can specify data once that all steppables respond to for consistent specification of units among all steppables
(rather than explicitly informing each steppable that you use in your simulation).

Currently, ``nCoVSteppableBase`` provides the following attributes that you can get or set on any instances of a class
that uses ``nCoVSteppableBase``, the changes to which will be reflected in all steppables that use
``nCoVSteppableBase``,

- ``step_period``: the period of one simulation step, in units of seconds per step (default 60).
- ``voxel_length``: the length of the side of each voxel, in units microns per voxel side length (default 1).

You can see an example of this in the default simulation specification distributed with the framework. The default
configuration specifies that the simulation consists of voxels of size 4 x 4 x 4 microns and five minutes per
simulation step. These operations can be found in ``Simulation/ViralInfectionVTM.py``,

.. code-block:: python

    from ViralInfectionVTMSteppables import CellInitializerSteppable
    steppable = CellInitializerSteppable(frequency=1)
    steppable.voxel_length = 4.0
    steppable.step_period = 5.0 * 60
    CompuCellSetup.register_steppable(steppable=steppable)

To be completely clear, you could perform these operations on any steppable that uses ``nCoVSteppableBase``, like
``CellInitializerSteppable`` in this example, and the changes will be reflected in all steppables so long as you make
them `before` calling ``CompuCelLSetup.run()``.
We will make every effort to ensure that all steppables in the shared library respond appropriately to such operations,
so that you don't need to worry about whether steppables are performing their model specification with the unit
specification that you prescribe.
We also recommend that you use these attributes in your own steppables that use ``nCoVSteppableBase``, and especially
should you decide to contribute your own model steppables to the shared library.

Framework Event System
----------------------

Say you've deployed a steppable from the shared library that creates a new cell, and say that you'd
like to do something with that cell when it is created. How would you know when it was created?
The computationally expensive solution is to check the cell inventory and look for changes.
However, some steppables in the shared library also issue procedures during such events, and tracking down where all
each steppable does such things can be tedious, confusing and even more computationally expensive.
Instead, ``nCoVSteppableBase`` provides an interface so that you can tell CC3D what additional things to do when such
events occur, whether they occur by a steppable from the shared library that you're using, or by one that you've
designed.
Likewise, all steppables in the shared library will respond to such events according to their model specification if
one of your steppables, or one of the shared library steppables you're using, issues them and is using
``nCoVSteppableBase``.

For each event, you can issue an event with a steppable use a particular function, and likewise your
steppable can respond to an event by defining a particular function, called a `callback`.
Currently, ``nCoVSteppableBase`` provides an interface to issue and respond to the events of

- creating a cell
- changing the type of a cell

Methods to issue these events and their callbacks are defined with the following interface,

.. code-block:: python

    class nCoVSteppableBase(SteppableBasePy):

        def new_cell(self, cell_type: int) -> CompuCell.CellG:
        '''
        Create and return a new cell
        `cell_type` is the integer id of the new cell's type according to a cc3d simulation.
        Equivalent to typical cc3d usage `new_cell(cell_type)`, with subsequent calls issued to `on_new_cell`
        for every registered `nCoVSteppableBase`-derived class by the framework.
        '''

        def on_new_cell(self, _new_cell: CompuCell.CellG) -> Union[None, bool]:
        '''
        A callback to respond to new_cell issued by a nCoVSteppableBase-derived class instance
        `_new_cell` is the newly created cell.
        '''

        def set_cell_type(self, cell: CompuCell.CellG, _type_id: int) -> None:
        '''
        Sets the type of a cell
        `cell` is the cell to which the change is made.
        `_type_id` is the integer id of the new type according to a cc3d simulation.
        Equivalent to typical cc3d usage `cell.type = _type_id`, with subsequent calls issued to `on_set_cell_type`
        for every registered `nCoVSteppableBase`-derived class by the framework.
        '''

        def on_set_cell_type(self, cell: CompuCell.CellG, old_type: int) -> Union[None, bool]:
        '''
        A callback to respond to set_cell_type issued by a nCoVSteppableBase-derived class instance.
        `cell` is the cell to which the change was made
        `old_type` is the previous type of the cell
        '''

For example, if one of your steppables uses ``nCoVSteppableBase`` and defines ``on_set_cell_type``, then
``on_set_cell_type`` will be called by the framework every time a steppable using ``nCoVSteppableBase`` changes the type
of a cell using ``set_cell_type``. Your steppable's implementation of ``on_set_cell_type`` can decide if the change in
cell type is relevant to your model specification, and if so, how to respond to it.
If your steppable isn't concerned with a particular event, then it can simply not define its corresponding callback.
Furthermore, all callbacks like ``on_new_cell`` and ``on_set_cell_type`` can also provide feedback to the framework
about if subsequent calls to the callback are needed.
For example, if your model specification requires information to decide about how to respond to an event that is not
yet available (`e.g.`, your steppable is waiting for information provided by another steppable's callback), your
steppable can notify the framework to call its callback again after calling the callback of every other registered
steppable that has not yet been called by returning ``False``.
If your callback does not require future calls, then it can return ``None``.
The ordering of calls to steppable callbacks is the same as the ordering of steppable registration with CC3D.

ODE Models in the Framework
---------------------------

``nCoVSteppableBase`` combines its event system with CC3D's built-in support for specifying, simulating and
manipulating ODE models defined in Antimony, CellML and SBML model syntax and attached to individual cells or defined
as free-floating (`i.e.`, simply running in the background).
``nCoVSteppableBase`` defines a method ``register_ode_model`` that registers an ODE model with the event system and
shares its information with all other registered steppables of a simulation that use ``nCoVSteppableBase``.
Likewise, any ODE model registered by a steppable from the shared libray will be available to your steppables if you
register the steppable from the shared library with CC3D, and the ODE model will also participate in the event system
and be simulated without any intervention by you or your steppables, but instead according to the specification of the
steppable from the shared library.

ODE models are registered with the framework as either free-floating, or as attached to particular cell type or set of
cell types.
When an ODE model is registered as free-floating, exactly one instance of the ODE model is automatically created and
shared with all registered steppables.
When an ODE model is registered as attached to a cell type or set of cell types, ODE model instances are automatically
instantiated, attached and destroyed by the framework during the creation of cells or changes to their type.
For example, if an ODE model is registered as corresponding to cell types "A" and "B", and a cell can have one of the
types "A", "B", "C" or "D", then the following events correspond to procedures performed by the framework,

- a cell of type "A" is created: the framework instantiates the ODE model and attaches it to the cell
- the cell changes to type "B": nothing occurs
- the cell changes to type "C": the ODE model is removed from the cell
- the cell changes to type "D": nothing occurs

Maintenace of the ODE models during an event is performed before issuing calls to event callbacks.

Typical CC3D usage of ODE models performs time integration of all ODE models with the method ``timestep_sbml``.
This functionality is strictly incompatible with this framework, since all ODE models are maintained by their parent
steppable, and so ``timestep_sbml`` should not be used unless no deployed model from the shared library registers an ODE
model. Otherwise, ``timestep_sbml`` issues time integration to all ODE models known by CC3D, which may result in
incorrect deployment of ODE models and simulation results.
Rather, steppables that use `nCoVSteppableBase` and register an ODE model with the framework should use the method
``timestep_ode_model`` to integrate the ODE model, which corresponds to calling ``timestep_sbml`` but for a single ODE
model.

All ODE models, whether attached to a cell or free-floating, can be accessed in using the typical CC3D fashion
(`e.g.`, ``self.sbml.MyODEModel``, ``cell.sbml.MyCellODEModel``, etc.).

The interface for registering, stepping and accesing ODE models is as follows,

.. code-block:: python

    class nCoVSteppableBase(SteppableBasePy):

        def register_ode_model(self,
                               model_name: str,
                               model_fcn: Callable,
                               ics_fcn: Callable = None,
                               cell_types: Union[str, Iterable[str], None] = None,
                               model_type: str = 'antimony',
                               step_size: float = 1.0) -> None:
        '''
        Registers an ode model with the framework and cc3d.
        `model_name` is the name of the ODE model
        `model_fcn` is a function that returns the string of the model when called;
            functions for free-floating models are passed no arguments;
            functions for cell-attached models are pass the cell to which the model is being attached
        `ics_fcn` is an optional function that returns a dictionary of initial conditions for the model
            functions for free-floating models are passed no arguments;
            functions for cell-attached models are pass the cell to which the model is being attached
        `cell_types` is an optional argument for specifying the corresponding cell type(s) associated with the model
        `model_type` specifies the language of the model; the default is antominy
        `step_size` is the time over which the model is integrated according to ODE model time in one integration step
        '''

        def timestep_ode_model(self, model_name: str) -> None:
        '''
        Integrate an ode model one step in time.
        `model_name` is the name of the model
        '''

        @property
        def ode_model_names(self) -> List[str]:
        '''
        List of ode models registered by all registered nCoVSteppableBase-based steppables
        '''

        def ode_models_by_cell_type(self, _cell_type: str):
        '''
        Returns a list of registered ode model names associated with a cell type name in a simulation.
        '''

        def cell_types_by_ode_model(self, _model_name: str) -> Union[None, List[int]]:
        '''
        Returns a list of cell type ids associated with a registered ode model in a simulation,
        or None if the model is free-floating
        '''

The framework adopts the convention that, for an event called by a method ``my_function``, there is a corresponding
callback ``on_my_function``.

Developing a Shared Module
==========================

Every module is defined with a unique name in the directory ``Models``. The space in the directory of your module is
your sandbox. There is no need to worry about colliding with developments by others, as it is your own unique space
within the greater framework.

When developing modules in this collection for usage in CC3D as add-on models, the main scripts defined in
the directory ``Simulation`` can be conveniently accessed in your module scripts (`e.g.`, for extending existing
models or accessing model inputs). The environment variable ``"ViralInfectionVTM"`` contains the path to the
root directory of the simulation framework. So, for example, to import the variable ``s_to_mcs`` from
``Simulation/ViralInfectionVTMModelInputs.py``, do the following basic Python procedures,

.. code-block:: python

    import os
    import sys
    sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
    from ViralInfectionVTMModelInputs import s_to_mcs

For a demonstration of this, see ``RecoverySteppables.py`` in the module ``Models.RecoverySimple``.

The same can be done for importing model modules defined in this collection (`e.g.`, for using or extending
add-on models). For example, if you would like to build a new steppable from ``SimpleRecoverySteppable`` defined
in the model module ``Models.RecoverySimple``, you can import ``SimpleRecoverySteppable`` with the following
basic Python procedures,

.. code-block:: python

    import os
    import sys
    sys.path.append(os.environ["ViralInfectionVTM"])
    from Models.RecoverySimple.RecoverySteppables import SimpleRecoverySteppable

Like any other Python class, steppables (and other code) defined in one model module can be extended by, or
integrated into, other modules. As such, the components of the overall simulation framework are not only
interchangable and shareable, but also `extensible`. For a demonstration of this, see ``RecoverySteppables.py``
in the module ``Models.RecoveryNeighbor``.

To promote shareability and extensibility, referenced CC3D data in a model specification should be implemented
dynamically with a clear API for how to tailor a steppable to a particular simulation.
For example, a model specification may be concerned with a particular cell type and field, each of which must be
assigned a name to be run in CC3D. If the cell type were named 'MyCellType" and the field were named 'MyField', then
typically CC3D specification would refer to each in a steppable using ``self.MYCELLTYPE`` and ``self.field.MyField``,
respectively.
Suppose that two steppables in two different modules describe different aspects of the same cell type and
field but so happen to name, and subsequently refer to, them differently.
This scenario would make it impossible for the user to use both modules and, hence, the two modules are incompatible.
As such, the internal CC3D data to which a steppable refers should be dynamically named, documented and configurable
through an API so that the user can inform module steppables of changes in cell type and field names.
The ability and methods to customize parameters and other internal data of a module is at the discretion of the module
developer.

For an example of dynamic naming, consider ``ViralInternalizationSteppable`` in
``Simulation.ViralInfectionVTMSteppables``, specifically concerning handling of the properties ``uninfected_type_name``
and ``target_field_name``, which define the names of the susceptible cell type and infecting name, respectively.
The steppable looks for a field named "Virus" by default, however if a user wanted to use a field named "InfluenzaA",
then they can do the following during import and registration of the steppable,

.. code-block:: python

    from ViralInfectionVTMSteppables import ViralInternalizationSteppable
    steppable = ViralInternalizationSteppable(frequency=1)
    steppable.target_field_name = "InfluenzaA"
    CompuCellSetup.register_steppable(steppable=steppable)

Module steppables should also define events and processes according to the aforementioned framework event system,
ODE model registration and framework-wide data.
A module that, for example, does not utilize ``set_cell_type`` will not be fully compatible with the overall framework,
and as such holds limited value in the shared library.
The same is true concerning framework-wide data like unit conversions, in that a module should incorporate all
information provided by the framework and adapt to user inputs appropriately according to the specification of the
module to maintain consistency with the simulation that user designs.
We will make every effort to provide feedback on how to modify code to accomplish this level of modularity in a
module, and also welcome comments and suggestions on how to better improve the framework and its documentation to
support easy incorporation of modules into the shared library.

Module Standards
----------------
Modules can be developed and incorporated into the framework shared library using standard GitHub practice of issuing
a pull request from a fork of the framework repository.
Pull requests will be reviewed according to the following standards.

All modules must maintain basic documentation. Documentation should be included as a multiline string in the module
``__init__.py`` with the following structure,

- A title heading for the module, followed by a basic overview of the module and any references to referred literature
- A heading "Maintainer(s)", with a list of all maintainers of the module and their affiliation(s)
- A heading "Contents", with a list of each file and directory and a brief description of its contents
- A heading "Change log", containing a sub-heading for each module version, each of which contains a list of changes

All modules must define and maintain in the module __init__.py the following current version information,

.. code-block:: python

    version_major = 0
    version_minor = 0
    version_build = 0
    version_str = f"{version_major}.{version_minor}.{version_build}"

Currently, versioning schemes are at the discretion of the module developer.

Each script must also provide at least a basic description of the contents of the script at the beginning of the
script, as appropriate for its type.

For CC3D-based model implementation, each steppable class definition must include a description of the steppable,
including its intended use, requirements and an overview of interacting with, manipulating and connecting it with
other modules of the framework.
Steppables should also direct the user with informative message when they are improperly used, so to help guide the
user on how to properly use the steppable in their simulation.
We also recommend, but do not enforce, a listing of every steppable class with a basic description at the beginning
of scripts that define multiple steppables.
For an example, see ``Simulation.ViralInfectionVTMSteppables.py``.

Currently no validation standards are enforced on modules, as the range of possible standards are too broad for a
community that includes both experimentalists and pure theoreticians.
However, we welcome module developers to refer to validation datasets and include subdirectories with scripts that
define validation routines as relevant to their module.

Support for additional functionality provided by the framework (`e.g.`, ``batchRun``) is at the discretion
of the module developer.
However, the framework development team reserves the right to incorporate modules into various functionality of the
framework.

Issues with, and suggestions for improvements to, the framework are welcome, and can be made as issues on the
repository of this framework.
