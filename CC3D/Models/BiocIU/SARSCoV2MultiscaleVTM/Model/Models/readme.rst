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
without overwriting each other (*e.g.*, ``Models.RecoverySimple.RecoverySteppables`` vs.
``Models.RecoveryNeighbor.RecoverySteppables``). The only necessarily unique aspect of a particular model
module is the name of its containing directory (*e.g.*, ``Models.RecoverySimple`` vs.
``Models.RecoveryNeighbor``). This scheme isolates model-specific development to the directory in which
the add-on model is defined, and modularizes the overall simulation framework into *shareable*,
*interchangable* model components.

When developing modules in this collection for usage in CC3D as add-on models, the main scripts defined in
the directory ``Simulation`` can be conveniently accessed in your module scripts (*e.g.*, for extending existing
models or accessing model inputs). The environment variable ``"ViralInfectionVTM"`` contains the path to the
root directory of the simulation framework. So, for example, to import the variable ``s_to_mcs`` from
``Simulation/ViralInfectionVTMModelInputs.py``, do the following basic Python procedures,

.. code-block:: python

    import os
    import sys
    sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
    from ViralInfectionVTMModelInputs import s_to_mcs

The same can be done for importing model modules defined in this collection (*e.g.*, for using or extending
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
interchangable and shareable, but also *extensible*.