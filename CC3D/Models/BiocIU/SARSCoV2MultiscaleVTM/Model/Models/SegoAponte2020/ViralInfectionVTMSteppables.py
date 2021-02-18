"""
Defines module steppables

Steppables
==========

CellsInitializerSteppable
-------------------------
Description: initializes epithelial sheet with a single infected cell

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import CellsInitializerSteppable

CompuCellSetup.register_steppable(steppable=CellsInitializerSteppable(frequency=1))

VirusFieldInitializerSteppable
------------------------------
Description: simple implementation of main VirusFieldInitializerSteppable to apply default virus field parameters

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import VirusFieldInitializerSteppable

CompuCellSetup.register_steppable(steppable=VirusFieldInitializerSteppable(frequency=1))

ViralReplicationSteppable
-------------------------
Description: implements viral replication module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ViralReplicationSteppable

CompuCellSetup.register_steppable(steppable=ViralReplicationSteppable(frequency=1))

ViralInternalizationSteppable
-----------------------------
Description: implements viral internalization module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ViralInternalizationSteppable

CompuCellSetup.register_steppable(steppable=ViralInternalizationSteppable(frequency=1))

ViralSecretionSteppable
-----------------------
Description: implements viral release module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ViralSecretionSteppable

CompuCellSetup.register_steppable(steppable=ViralSecretionSteppable(frequency=1))

ImmuneCellKillingSteppable
--------------------------
Description: implements immune cell direct cytotoxicity and bystander effect modules

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ImmuneCellKillingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellKillingSteppable(frequency=1))

ChemotaxisSteppable
-------------------
Description: implements immune cell chemotaxis module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ChemotaxisSteppable

CompuCellSetup.register_steppable(steppable=ChemotaxisSteppable(frequency=1))

ImmuneCellSeedingSteppable
--------------------------
Description: implements immune cell seeding and removal of immune cell recruitment module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ImmuneCellSeedingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellSeedingSteppable(frequency=1))

SimDataSteppable
----------------
Description: plots/writes simulation data of interest

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import SimDataSteppable

CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))

CytokineProductionAbsorptionSteppable
-------------------------------------
Description: implements cytokine production/secretion and immune cell activation module, and sets cytokine field data

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import CytokineProductionAbsorptionSteppable

CompuCellSetup.register_steppable(steppable=CytokineProductionAbsorptionSteppable(frequency=1))

ImmuneRecruitmentSteppable
--------------------------
Description: implements immune cell recruitment module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ImmuneRecruitmentSteppable

CompuCellSetup.register_steppable(steppable=ImmuneRecruitmentSteppable(frequency=1))

oxidationAgentModelSteppable
----------------------------
Description: implements immune cell oxidizing agent cytotoxicity module

Usage: In ViralInfectionVTM.py, add the following

from Models.SegoAponte2020.ViralInfectionVTMSteppables import oxidationAgentModelSteppable

CompuCellSetup.register_steppable(steppable=oxidationAgentModelSteppable(frequency=1))
"""

import math
import numpy as np

from . import ViralInfectionVTMBasePy
from . import ViralInfectionVTMLib
from . import ViralInfectionVTMModelInputs

# Import main steppables module
import ViralInfectionVTMSteppables as MainSteppables

# Import toolkit
from nCoVToolkit import nCoVUtils

ViralInfectionVTMSteppableBasePy = ViralInfectionVTMBasePy.ViralInfectionVTMSteppableBasePy

# Name of generic immune cell type
immune_type_name = 'Immunecell'
# Name of generic cytokine field
cytokine_field_name = 'cytokine'
# Name of generic oxidative agent field
oxidator_field_name = 'oxidator'

rng = np.random  # alias for random number generators (rng)


class CellsInitializerSteppable(ViralInfectionVTMSteppableBasePy, MainSteppables.CellInitializerSteppable):
    """
    Initializes epithelial sheet with a single infected cell

    Requires ViralReplicationSteppable registered beforehand to automatically instantiate cells with viral
    replication model. Otherwise, viral replication model needs to be instantiated per cell after initialization
    """

    unique_key = ViralInfectionVTMLib.cell_initializer_steppable_key

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        MainSteppables.CellInitializerSteppable.__init__(self, frequency)

        # Initialize default data
        self.set_cell_diameter(ViralInfectionVTMModelInputs.cell_diameter)
        self.set_volume_parameter(ViralInfectionVTMModelInputs.volume_lm)
        self.single_infected_cell()

    def start(self):
        """
        Called once to initialize simulation
        """
        ViralInfectionVTMSteppableBasePy.start(self)
        MainSteppables.CellInitializerSteppable.start(self)

        for cell in self.cell_list_by_type(self.uninfected_type_id, self.infected_type_id):
            if ViralInfectionVTMLib.vrl_key in cell.dict.keys():
                ViralInfectionVTMLib.pack_viral_replication_variables(cell)

        # Set initially infected cell data
        for cell in self.cell_list_by_type(self.infected_type_id):
            var_unpacking = ViralInfectionVTMLib.vr_cell_dict_to_sym[ViralInfectionVTMLib.vrm_unpacking]
            getattr(cell.sbml, self.vr_model_name)[var_unpacking] = 1.0


class VirusFieldInitializerSteppable(MainSteppables.VirusFieldInitializerSteppable):
    """Simple implementation of main VirusFieldInitializerSteppable to apply default virus field parameters"""

    unique_key = ViralInfectionVTMLib.virus_field_initializer_key

    def __init__(self, frequency=1):
        super().__init__(frequency=frequency)

        # Initialize default data
        self.set_diffusion_coefficient(ViralInfectionVTMModelInputs.virus_dc)
        self.set_decay_coefficient(ViralInfectionVTMModelInputs.virus_decay)


class ViralReplicationSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral replication module

    Assigns a callback to instantiation of new epithelial cells to assign a new intracellular model

    Additional cell types can be registered and unregistered with `register_type` and `unregister_type`, respectively.
    Cell types that are registered with the keyword argument `infected` set to True will have their intracellular model
    integrated in time every simulation step. All registrations must occur before start is called.

    By default, module uninfected, infected and virus-releasing types are registered, with the latter two registered as
    infected.

    The criterion for the transition from module types infected to virus releasing and virus releasing to dead are
    described by `cell_releases` and `cell_dies`, respectively, which can be set/overridden as desired.

    The step size of the viral replication model per simulation step can be set with the attribute `step_size`.
    This is used as a reference value for calculations, and does not set the step size used in the actual viral
    replication model.

    The threshold at which cells convert from the infected to virus-releasing type can be set with the attribute
    `cell_infection_threshold`.

    The dissociation and Hill coefficients for cell death can be set with the attributes `diss_coeff_death` and
    `hill_coeff_death`, respectively.

    Model parameters for new viral replication model instantiations can be set with `set_vrm_param`. Existing models
    must be modified by referencing the instance of the existing model in the typical way.

    Reports the occurrence of death if the attribute `simdata_steppable` is set or a module `SimDataSteppable` instance
    is found in the shared dictionary.
    If `simdata_steppable` is set, then the referred instance must have the method `track_death_viral`, which it will
    call without arguments when reporting death.
    """

    unique_key = ViralInfectionVTMLib.vrm_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self._registered_types = []
        self._infected_types = []
        self._step_size = 1.0
        self._cell_infection_threshold = 1.0
        self._diss_coeff_uptake_apo = 1.0
        self._hill_coeff_uptake_apo = 0.0
        self._vrm_params = {'unpacking_rate': 0.0,
                            'replicating_rate': 0.0,
                            'r_half': 0.0,
                            'translating_rate': 0.0,
                            'packing_rate': 0.0,
                            'secretion_rate': 0.0}

        # Initialize default data
        self.sbml_options = {'relative': 1e-10, 'absolute': 1e-12}
        self.set_step_size(ViralInfectionVTMModelInputs.vr_step_size)
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name, infected=True)
        self.register_type(MainSteppables.virus_releasing_type_name, infected=True)
        self.set_cell_infection_threshold(ViralInfectionVTMModelInputs.cell_infection_threshold)
        self.set_diss_coeff_death(ViralInfectionVTMModelInputs.diss_coeff_uptake_apo)
        self.set_hill_coeff_death(ViralInfectionVTMModelInputs.hill_coeff_uptake_apo)
        self.set_vrm_param(unpacking_rate=ViralInfectionVTMModelInputs.unpacking_rate,
                           replicating_rate=ViralInfectionVTMModelInputs.replicating_rate,
                           r_half=ViralInfectionVTMModelInputs.r_half,
                           translating_rate=ViralInfectionVTMModelInputs.translating_rate,
                           packing_rate=ViralInfectionVTMModelInputs.packing_rate,
                           secretion_rate=ViralInfectionVTMModelInputs.secretion_rate)

    def start(self):
        """
        Called once to initialize simulation
        """
        # Load model
        self.set_sbml_global_options(self.sbml_options)

        self.register_ode_model(model_name=self.vr_model_name,
                                model_fcn=self.epithelial_model_fcn_generator(),
                                cell_types=self._registered_types,
                                step_size=self._step_size)

    def step(self, mcs):
        """
        Called every simulation step
        """
        if self.simdata_steppable is None:
            try:
                self.simdata_steppable: SimDataSteppable = self.shared_steppable_vars[SimDataSteppable.unique_key]
            except KeyError:
                pass

        # Do viral model
        for cell in self.cell_list_by_type(*self.infected_type_ids):
            # Step the model for this cell
            ViralInfectionVTMLib.step_sbml_model_cell(cell=cell)
            # Pack state variables into cell dictionary
            ViralInfectionVTMLib.pack_viral_replication_variables(cell=cell)

            # Test for infected -> virus-releasing
            if self.cell_releases(cell):
                self.set_cell_type(cell, self.virus_releasing_type_id)

            # Test for virus-releasing -> dead
            if self.cell_dies(cell):
                self.kill_cell(cell=cell)
                if self.simdata_steppable is not None:
                    self.simdata_steppable.track_death_viral()

    def cell_releases(self, _cell):
        """
        Criterion describing transition to virus-releasing

        :param _cell: a cell
        :return: True if cell releases
        """
        return _cell.type == self.infected_type_id and \
               _cell.dict[ViralInfectionVTMLib.vrm_assembled] > self._cell_infection_threshold

    def cell_dies(self, _cell):
        """
        Criterion describing transition to dead

        :param _cell: a cell
        :return: True if cell dies
        """
        return _cell.type == self.virus_releasing_type_id and \
               np.random.random() < nCoVUtils.hill_equation(_cell.dict[ViralInfectionVTMLib.vrm_assembled],
                                                            self._diss_coeff_uptake_apo,
                                                            self._hill_coeff_uptake_apo)

    def epithelial_model_fcn_generator(self):
        """
        Get model string generator callback
        """
        def model_fcn(_cell):
            """
            Callback to the model string generator
            """
            _cell.dict[ViralInfectionVTMLib.vrl_key] = True
            for k in [ViralInfectionVTMLib.vrm_unpacking,
                      ViralInfectionVTMLib.vrm_replicating,
                      ViralInfectionVTMLib.vrm_packing,
                      ViralInfectionVTMLib.vrm_assembled,
                      ViralInfectionVTMLib.vrm_uptake]:
                if k not in _cell.dict.keys():
                    _cell.dict[k] = 0.0
            return self.viral_replication_model_string(self._vrm_params['unpacking_rate'] * self.step_period,
                                                       self._vrm_params['replicating_rate'] * self.step_period,
                                                       self._vrm_params['r_half'],
                                                       self._vrm_params['translating_rate'] * self.step_period,
                                                       self._vrm_params['packing_rate'] * self.step_period,
                                                       self._vrm_params['secretion_rate'] * self.step_period,
                                                       _cell.dict[ViralInfectionVTMLib.vrm_unpacking],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_replicating],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_packing],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_assembled],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_uptake])
        return model_fcn

    def register_type(self, _type_name: str, infected: bool = False):
        """
        Register a cell type for simulation of viral replication model

        :param _type_name: name of the cell type
        :param infected: flag describing if cell releases virus
        :return: None
        """
        self._registered_types.append(_type_name)
        if infected:
            self._infected_types.append(_type_name)

    def unregister_type(self, _type_name: str):
        """
        Unregister a cell type for simulation of viral replication model

        :param _type_name: name of cell type
        :return: None
        """
        self._registered_types.remove(_type_name)
        if _type_name in self._infected_types:
            self._infected_types.remove(_type_name)

    @property
    def infected_type_ids(self):
        return [getattr(self, x.upper()) for x in self._infected_types]

    @property
    def registered_type_ids(self):
        return [getattr(self, x.upper()) for x in self._registered_types]

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback. Ensures that necessary data is available and properly set
        """
        if _new_cell.type in self.registered_type_ids and ViralInfectionVTMLib.vrl_key not in _new_cell.dict.keys():
            _new_cell.dict[ViralInfectionVTMLib.vrl_key] = False

    def set_step_size(self, _val: float):
        """
        Set viral replication model step size

        :param _val: step size
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Step size must be positive")
        self._step_size = _val

    @property
    def step_size(self):
        """
        Threshold above which cells become virus-releasing
        """
        return self._step_size

    @step_size.setter
    def step_size(self, _val: float):
        self.set_step_size(_val)

    def set_cell_infection_threshold(self, _val: float):
        """
        Set threshold above which cells become virus-releasing

        :param _val: infection threshold
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Cell infection threshold must be positive")
        self._cell_infection_threshold = _val

    @property
    def cell_infection_threshold(self):
        return self._cell_infection_threshold

    @cell_infection_threshold.setter
    def cell_infection_threshold(self, _val: float):
        self.set_cell_infection_threshold(_val)

    def set_diss_coeff_death(self, _val: float):
        """
        Set dissociation coefficient for virally-induced death

        :param _val: dissociation coefficient
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Dissociation coefficient must be positive")
        self._diss_coeff_uptake_apo = _val

    @property
    def diss_coeff_death(self):
        """
        Dissociation coefficient for virally-induced death
        """
        return self._diss_coeff_uptake_apo

    @diss_coeff_death.setter
    def diss_coeff_death(self, _val):
        self.set_diss_coeff_death(_val)

    def set_hill_coeff_death(self, _val: float):
        """
        Set Hill coefficient for virally-induced death

        :param _val: Hill coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Hill coefficient must be non-negative")
        self._hill_coeff_uptake_apo = _val

    @property
    def hill_coeff_death(self):
        """
        Hill coefficient for virally-induced death
        """
        return self._hill_coeff_uptake_apo

    @hill_coeff_death.setter
    def hill_coeff_death(self, _val):
        self.set_hill_coeff_death(_val)

    def set_vrm_param(self, **kwargs):
        """
        Set parameters in viral replication model for newly created cells

        Keywords are

        - unpacking_rate (units 1/s)
        - replicating_rate (units 1/s)
        - r_half (ul)
        - translating_rate (units 1/s)
        - packing_rate (units 1/s)
        - secretion_rate (units 1/s)

        Any model parameter can be set, so long as the value is non-negative

        :param kwargs: keyword argument values
        :return: None
        """
        for k, v in kwargs.items():
            if k not in self._vrm_params.keys():
                raise AttributeError(f'Unrecognized parameter ({k} = {v})')
            elif v < 0.0:
                raise ValueError(f'Value must be non-negative ({k} = {v})')
            self._vrm_params[k] = v

    def get_vrm_param(self) -> dict:
        """
        Get a copy of parameters in viral replication model for newly created cells
        """
        return self._vrm_params.copy()


class ViralInternalizationSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral internalization module

    All new cells of the uninfected type are initialized with a number of unbound receptors, which is stored in the
    cell dictionary with a key of the value of `unbound_receptors_cellg_key` defined in `ViralInfectionVTMLib`.
    Updates to the number of unbound receptors on a cell's surface can be made with `update_cell_receptors`.

    The number of unbound receptors applied during initialization can be set with `initial_unbound_receptors`. Note that
    this value is also used during internalization calculations.

    When internalization occurs, an uninfected cell's type is changed to infected.
    When the number of unbound receptors is zero, internalization is not considered.

    Then name of the uninfected and infected types can be set with `uninfected_type_name` and `infected_type_name`,
    respectively.

    Model parameters used in the calculations made in the call to `do_cell_internalization` can be set with attributes
    as follows

    - Number of initial unbound receptors on newly created cells: `initial_unbound_receptors`
    - Virus-receptor association affinity: `kon`
    - Virus-receptor dissociation affinity: `koff`
    - Internalization Hill coefficient: `hill_coeff_uptake`
    - Uptake rate coefficient: `rate_coeff_uptake`

    """

    unique_key = ViralInfectionVTMLib.vim_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self._initial_unbound_receptors = 0.0
        self._kon = 1.0
        self._koff = 1.0
        self._hill_coeff_uptake = 1.0
        self._rate_coeff_uptake = 1.0

        # Initialize default data
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_initial_unbound_receptors(ViralInfectionVTMModelInputs.initial_unbound_receptors)
        self.set_kon(ViralInfectionVTMModelInputs.kon)
        self.set_koff(ViralInfectionVTMModelInputs.koff)
        self.set_hill_coeff_uptake(ViralInfectionVTMModelInputs.hill_coeff_uptake_pr)
        self.set_rate_coeff_uptake(ViralInfectionVTMModelInputs.rate_coeff_uptake_pr)

    def do_cell_internalization(self, cell, viral_amount_com):
        """
        Describes internalization criterion

        :param cell: a cell
        :param viral_amount_com: total amount of virus associated with a cell
        :return: boolean of whether internalization occurs, and corresponding uptake amount
        :rtype: (bool, float)
        """
        if cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] == 0:
            return False, 0.0

        # Conversion to cc3d units
        pmol_to_cc3d_au = ViralInfectionVTMModelInputs.pmol_to_cc3d_au
        kon = self._kon * self.step_period / self.voxel_length**3 / pmol_to_cc3d_au
        koff = self._koff * self.step_period

        _k = kon * cell.volume / koff
        hill_coeff_uptake_pr = self._hill_coeff_uptake
        receptors = cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key]
        diss_coeff_uptake_pr = (self._initial_unbound_receptors / 2.0 / _k / receptors) ** (1.0 / hill_coeff_uptake_pr)
        uptake_probability = nCoVUtils.hill_equation(viral_amount_com,
                                                     diss_coeff_uptake_pr,
                                                     hill_coeff_uptake_pr)

        cell_does_uptake = np.random.rand() < uptake_probability
        uptake_amount = self.step_period / self._rate_coeff_uptake * uptake_probability

        if cell_does_uptake and cell.type == self.uninfected_type_id:
            self.set_cell_type(cell, self.infected_type_id)

        return cell_does_uptake, uptake_amount

    def update_cell_receptors(self, cell, receptors_increment):
        """
        Updates the number of unbound receptors on a cell's surface

        :param cell: a cell
        :param receptors_increment: increment to number of receptors
        :return: None
        """
        cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] = max(
            cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] + receptors_increment, 0.0)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback. Ensures that necessary data is available and properly set
        """
        if _new_cell.type == self.uninfected_type_id and \
                ViralInfectionVTMLib.unbound_receptors_cellg_key not in _new_cell.dict.keys():
            _new_cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] = self._initial_unbound_receptors

    def set_initial_unbound_receptors(self, _val: float):
        """
        Set number of initial unbound receptors on newly created cells

        :param _val: number of initial unbound receptors
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Initial unbound receptors must be positive")
        self._initial_unbound_receptors = _val

    @property
    def initial_unbound_receptors(self):
        """
        Number of initial unbound receptors
        """
        return self._initial_unbound_receptors

    @initial_unbound_receptors.setter
    def initial_unbound_receptors(self, _val: float):
        self.set_initial_unbound_receptors(_val)

    def set_kon(self, _val: float):
        """
        Set virus-receptor association affinity, in units 1/(M * s)

        :param _val: Virus-receptor association affinity
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("kon must be positive")
        self._kon = _val

    @property
    def kon(self):
        """
        Virus-receptor association affinity, in units 1/(M * s)
        """
        return self._kon

    @kon.setter
    def kon(self, _val: float):
        self.set_kon(_val)

    def set_koff(self, _val: float):
        """
        Set virus-receptor dissociation affinity, in units 1/s

        :param _val: virus-receptor dissociation affinity
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("koff must be positive")
        self._koff = _val

    @property
    def koff(self):
        """
        Virus-receptor dissociation affinity, in units 1/s
        """
        return self._koff

    @koff.setter
    def koff(self, _val: float):
        self.set_koff(_val)

    def set_hill_coeff_uptake(self, _val: float):
        """
        Set internalization Hill coefficient

        :param _val: internalization Hill coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Hill coefficient must non-negative")
        self._hill_coeff_uptake = _val

    @property
    def hill_coeff_uptake(self):
        """
        Internalization Hill coefficient
        """
        return self._hill_coeff_uptake

    @hill_coeff_uptake.setter
    def hill_coeff_uptake(self, _val: float):
        self.set_hill_coeff_uptake(_val)

    def set_rate_coeff_uptake(self, _val):
        """
        Set uptake rate coefficient, in units s

        :param _val: uptake rate coefficient
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Uptake rate coefficient must be positive")
        self._rate_coeff_uptake = _val

    @property
    def rate_coeff_uptake(self):
        """
        Uptake rate coefficient, in units s
        """
        return self._rate_coeff_uptake

    @rate_coeff_uptake.setter
    def rate_coeff_uptake(self, _val: float):
        self.set_rate_coeff_uptake(_val)


class ViralSecretionSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral release module

    This steppable uses a callback per registered cell type that is called on each cell of the registered type every
    time step to calculate interactions with the module virus field.

    Every secretion callback has the signature `secr_func(_steppable, _cell, _mcs)`

    - `_steppable` (ViralSecretionSteppable): self
    - `_cell` (cc3d.cpp.CompuCell.CellG): a cell
    - `_mcs` (int): the current simulation step

    Callbacks can be registered for a cell type with `register_secretion_by_type`, and unregistered for a cell type with
    `unregister_secretion_by_type`.

    Requires a virus-releasing cell type, the name of which can be set with the attribute `virus_releasing_type_name`.

    Requires the module viral replication module as managed by `ViralReplicationSteppable`,
    or one of the same name that defines the variable `secretion_rate`.

    The rate of viral release assigned to each a virus-releasing cell when changing type can be set with the attribute
    `viral_secretion_rate`.

    By default, module uninfected, infected and virus-releasing types are registered with the callback
    `secr_func_epithelial`.

    By default, requires a diffusion field with module virus field name. If all types registered with the callback
    `secr_func_epithelial` are unregistered, then this requirement is no longer relevant.
    The name of the virus field can be set with the attribute `virus_field_name`.

    By default, requires `ViralInternalizationSteppable` for determining/implementing internalization events. This can
    be customized by setting `vim_steppable` with any object that implements `do_cell_internalization` and
    `update_cell_receptors` with the same signature as those defined by `ViralInternalizationSteppable`. If all types
    registered with the callback `secr_func_epithelial` are unregistered, then this requirement is no longer relevant.
    """

    unique_key = ViralInfectionVTMLib.vrs_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        # Reference to ViralInternalizationSteppable
        self.vim_steppable = None

        self._secretors_by_type = {}
        self._secretion_rate = 0.0

        # Initialize default data
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_virus_field_name(MainSteppables.virus_field_name)
        self.register_secretion_by_type(MainSteppables.uninfected_type_name, self.secr_func_epithelial)
        self.register_secretion_by_type(MainSteppables.infected_type_name, self.secr_func_epithelial)
        self.register_secretion_by_type(MainSteppables.virus_releasing_type_name, self.secr_func_epithelial)
        self.set_viral_secretion_rate(ViralInfectionVTMModelInputs.secretion_rate)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        if self.vim_steppable is None:
            try:
                self.vim_steppable: ViralInternalizationSteppable = \
                    self.shared_steppable_vars[ViralInternalizationSteppable.unique_key]
            except KeyError:
                pass

        for type_name, secr_func in self._secretors_by_type.items():
            for cell in self.cell_list_by_type(getattr(self, type_name.upper())):
                secr_func(self, cell, mcs)

    @staticmethod
    def secr_func_epithelial(self, _cell, _mcs):
        """
        Secretion function callback for cells registered as epithelial types by default
        """
        # Evaluate probability of cell uptake of viral particles from environment
        # If cell isn't infected, it changes type to infected here if uptake occurs
        secretor = self.get_field_secretor(self._virus_field_name)
        viral_amount_com = self.total_seen_field(field=self.virus_field, cell=_cell, estimate=True)
        cell_does_uptake, uptake_amount = self.vim_steppable.do_cell_internalization(_cell, viral_amount_com)
        if cell_does_uptake:
            uptake = secretor.uptakeInsideCellTotalCount(_cell, 1E12, uptake_amount / _cell.volume)
            _cell.dict[ViralInfectionVTMLib.vrm_uptake] = abs(uptake.tot_amount)
            self.vim_steppable.update_cell_receptors(
                cell=_cell,
                receptors_increment=-_cell.dict[ViralInfectionVTMLib.vrm_uptake]*self.step_period)
            ViralInfectionVTMLib.set_viral_replication_cell_uptake(cell=_cell,
                                                                   uptake=_cell.dict[ViralInfectionVTMLib.vrm_uptake])

        if _cell.type == self.virus_releasing_type_id:
            sec_amount = ViralInfectionVTMLib.get_viral_replication_cell_secretion(cell=_cell)
            secretor.secreteInsideCellTotalCount(_cell, sec_amount / _cell.volume)

    def register_secretion_by_type(self, _type_name: str, _secr_func):
        """
        Register a secretion callback by cell type

        :param _type_name: name of cell type
        :param _secr_func: secretion callback
        :return: None
        """
        self._secretors_by_type[_type_name] = _secr_func

    def unregister_secretion_by_type(self, _type_name: str):
        """
        Unregister a secretion callback by cell type

        :param _type_name: name of cell type
        :return: callback of registered cell type
        """
        return self._secretors_by_type.pop(_type_name)

    def set_viral_secretion_rate(self, _val: float):
        """
        Set viral release rate for virus-releasing cell type, in units 1/s

        :param _val: viral release rate
        :return: None
        """
        self._secretion_rate = _val

    @property
    def viral_secretion_rate(self):
        """
        Viral release rate for virus-releasing cell type, in units 1/s
        """
        return self._secretion_rate

    @viral_secretion_rate.setter
    def viral_secretion_rate(self, _val: float):
        self.set_viral_secretion_rate(_val)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback
        """
        self.on_set_cell_type(cell=_new_cell, old_type=None)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback. If a cell's type changes to the registered virus-releasing type, then
        release of virus in enabled
        """
        if cell.type == self.virus_releasing_type_id:
            ViralInfectionVTMLib.enable_viral_secretion(cell=cell,
                                                        secretion_rate=self._secretion_rate * self.step_period)


class ImmuneCellKillingSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell direct cytotoxicity and bystander effect modules.

    The name of the immune cell type can be set with the attribute `immune_type_name`.

    Cell types that can be killed by contact-mediated interactions with the module immune cell type can be
    registered and unregistered with `append_contact_target` and `remove_contact_target`, respectively.

    Cell types that can be killed by the bystander effect can be registered and unregistered with `add_bystander_target`
    and `remove_bystander_target`, respectively.

    By default, module infected and virus-releasing types are registered for contact-mediated killing.

    By default, module uninfected, infected and virus-releasing types are registered for bystander effect killing.

    Reports the occurrence of death by contact and bystander effect if the attribute `simdata_steppable` is set or
    a module `SimDataSteppable` instance is found in the shared dictionary.
    If `simdata_steppable` is set, then the referred instance must have the methods `track_death_contact` and
    `track_death_bystander`, which it will call without arguments when reporting contact and bystander effect death,
    respectively.
    """

    unique_key = ViralInfectionVTMLib.immune_killing_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self._contact_targets = []
        self._bystander_targets = {}

        # Initialize defaults
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.set_immune_type_name(immune_type_name)
        # Infected and virus-releasing types are contact targets
        self.append_contact_target(MainSteppables.infected_type_name)
        self.append_contact_target(MainSteppables.virus_releasing_type_name)
        # Uninfected, infected and virus-releasing types are bystander targets
        self.add_bystander_target(MainSteppables.uninfected_type_name, ViralInfectionVTMModelInputs.bystander_effect)
        self.add_bystander_target(MainSteppables.infected_type_name, ViralInfectionVTMModelInputs.bystander_effect)
        self.add_bystander_target(MainSteppables.virus_releasing_type_name,
                                  ViralInfectionVTMModelInputs.bystander_effect)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        if self.simdata_steppable is None:
            try:
                self.simdata_steppable: SimDataSteppable = self.shared_steppable_vars[SimDataSteppable.unique_key]
            except KeyError:
                pass

        killed_cells = []
        for cell in self.cell_list_by_type(*self.contact_target_ids):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.immune_type_id:
                        self.kill_cell(cell=cell)
                        killed_cells.append(cell)
                        if self.simdata_steppable is not None:
                            self.simdata_steppable.track_death_contact()

        # Bystander Effect
        for cell in killed_cells:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type in self.bystander_target_ids:
                        if np.random.random() < self._bystander_targets[self.get_type_name_by_cell(neighbor)]:
                            self.kill_cell(cell=neighbor)
                            if self.simdata_steppable is not None:
                                self.simdata_steppable.track_death_bystander()

    def append_contact_target(self, _name: str):
        """
        Add a cell type target for contact-mediated killing

        :param _name: name of cell type
        :return: None
        """
        self._contact_targets.append(_name)

    def remove_contact_target(self, _name: str):
        """
        Remove a cell type target for contact-mediated killing

        :param _name: name of cell type
        :return: None
        """
        self._contact_targets.remove(_name)

    def contact_target_names(self) -> list:
        """
        Names of contact targets
        """
        return self._contact_targets.copy()

    @property
    def contact_target_ids(self):
        """
        Ids of contact targets according to a cc3d simulation
        """
        return [getattr(self, x.upper()) for x in self.contact_target_names()]

    def add_bystander_target(self, _name: str, _prob: float):
        """
        Add a cell type target for bystander effect killing

        :param _name: name of cell type
        :param _prob: probability of killing by the bystander effect
        :return: None
        """
        if not 0 <= _prob <= 1:
            raise ValueError('Bystander effect probability must be in [0, 1]')
        self._bystander_targets[_name] = _prob

    def remove_bystander_target(self, _name: str):
        """
        Remove a cell type target for bystander effect killing

        :param _name: name of cell type
        :return: None
        """
        self._bystander_targets.pop(_name)

    @property
    def bystander_targets(self):
        """
        Names of bystander effect targets
        """
        return [x for x in self._bystander_targets.keys()]

    @property
    def bystander_target_ids(self):
        """
        Ids of bystander effect targets according to a cc3d simulation
        """
        return [getattr(self, x.upper()) for x in self.bystander_targets]

    @property
    def bystander_target_data(self):
        """
        Bystander effect data

        Keys are names of register cell types, items are probabilities of death
        """
        return self._bystander_targets.copy()


class ChemotaxisSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell chemotaxis module

    By default, immune cells chemotaxis according to the module cytokine field.
    The name of the field can be set with the attribute `target_field_name`.

    The name of the immune cell type can be set with the attribute `immune_type_name`.

    The chemotaxis parameter of immune cell types can be set with the attribute `lamda_chemotaxis`.
    """

    unique_key = ViralInfectionVTMLib.chemotaxis_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self._target_field_name = ''
        self._lamda_chemotaxis = 0.0

        # Initialize defaults
        self.set_immune_type_name(immune_type_name)
        self.set_target_field_name(cytokine_field_name)
        self.set_lamda_chemotaxis(ViralInfectionVTMModelInputs.lamda_chemotaxis)

    def start(self):
        """
        Called once to initialize simulation
        """
        for cell in self.cell_list_by_type(self.immune_type_id):
            self.add_cell_chemotaxis_data(cell=cell)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        field = self.target_field
        for cell in self.cell_list_by_type(self.immune_type_id):

            cd = self.chemotaxisPlugin.getChemotaxisData(cell, self._target_field_name)
            if cell.dict[ViralInfectionVTMLib.activated_cellg_key]:
                cd.setLambda(self._lamda_chemotaxis / (1.0 + field[cell.xCOM, cell.yCOM, 1]))
            else:
                cd.setLambda(0)

    def set_target_field_name(self, _name: str):
        """
        Set name of chemotaxis target field

        :param _name: name of chemotaxis target field
        :return: None
        """
        self._target_field_name = _name

    @property
    def target_field_name(self):
        """
        Name of chemotaxis target field
        """
        return self._target_field_name

    @target_field_name.setter
    def target_field_name(self, _name: str):
        self.set_target_field_name(_name)

    def set_lamda_chemotaxis(self, _val: float):
        """
        Set chemotaxis parameter value

        :param _val: chemotaxis parameter value
        :return: None
        """
        self._lamda_chemotaxis = _val

    @property
    def lamda_chemotaxis(self):
        """
        Chemotaxis parameter value
        """
        return self._lamda_chemotaxis

    @lamda_chemotaxis.setter
    def lamda_chemotaxis(self, _val: float):
        self.set_lamda_chemotaxis(_val)

    @property
    def target_field(self):
        """
        Reference to chemotaxis target field
        """
        return getattr(self.field, self._target_field_name)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback.

        If a cell is of registered immune cell type, then its chemotaxis data is initialized
        """
        if _new_cell.type == self.immune_type_id:
            self.add_cell_chemotaxis_data(cell=_new_cell)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback.

        If a cell is of registered immune cell type, then its chemotaxis data is initialized
        """
        if cell.type == self.immune_type_id:
            self.add_cell_chemotaxis_data(cell=cell)

    def add_cell_chemotaxis_data(self, cell):
        """
        Add chemotaxis data to a cell

        :param cell: a cell
        :return: None
        """
        cd = self.chemotaxisPlugin.addChemotaxisData(cell, self._target_field_name)
        cd.assignChemotactTowardsVectorTypes([self.MEDIUM])


class ImmuneCellSeedingSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell seeding and removal of immune cell recruitment module.

    Seeded immune cells are initialized with the model parameters of this module in their dictionary.

    The name of the field by which immune cells are seeded can be set with the attribute `seeding_field_name`.

    The target volume and volume parameter of newly created immune cells can be set with the attributes
    `cell_volume` and `volume_parameter`, respectively.

    By default, uses `ImmuneRecruitmentSteppable` to determine seeding and removal probabilities.
    If customizing, this functionality can be replaced by setting the attribute `ir_steppable` with any object that
    provides methods `get_immune_seeding_prob` and `get_immune_removal_prob` that return the probability of seeding and
    removing an immune cell, respectively.
    If deriving, this functionality can be replaced by overriding `get_immune_seeding_prob` and
    `get_immune_removal_prob` to describe the tests for seeding and removing an immune cell, respectively.
    """

    unique_key = ViralInfectionVTMLib.immune_seeding_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to ImmuneResponseSteppable
        self.ir_steppable = None

        self._seeding_field_name = ''
        self._cell_volume = 0
        self._volume_lm = 0
        self._initial_immune_seeding = 0

        # Initialize defaults
        self.set_immune_type_name(immune_type_name)
        self.set_seeding_field_name(MainSteppables.virus_field_name)
        self.set_cell_volume(ViralInfectionVTMModelInputs.cell_volume)
        self.set_volume_parameter(ViralInfectionVTMModelInputs.volume_lm)
        self.set_initial_immune_seeding(ViralInfectionVTMModelInputs.initial_immune_seeding)

    def start(self):
        """
        Called once to initialize simulation
        """

        cell_diameter = int(self.cell_diameter / self.voxel_length)
        for iteration in range(int(self._initial_immune_seeding)):
            cell = True
            while cell:
                xi = np.random.randint(0, self.dim.x - 2 * cell_diameter)
                yi = np.random.randint(0, self.dim.y - 2 * cell_diameter)
                for x in range(xi, xi + cell_diameter):
                    for y in range(yi, yi + cell_diameter):
                        cell = self.cell_field[x, y, 1]
                        break
                cell = False
            cell = self.new_cell(self.immune_type_id)
            self.cell_field[x:x + cell_diameter, y:y + cell_diameter, 1] = cell

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        if self.ir_steppable is None:
            try:
                self.ir_steppable: ImmuneRecruitmentSteppable = \
                    self.shared_steppable_vars[ImmuneRecruitmentSteppable.unique_key]
            except KeyError:
                pass

        for cell in self.cell_list_by_type(self.immune_type_id):
            p_immune_dying = np.random.random()
            if p_immune_dying < self.get_immune_removal_prob():
                cell.targetVolume = 0.0

        p_immune_seeding = np.random.random()
        if p_immune_seeding < self.get_immune_seeding_prob():
            open_space = True
            viral_concentration = 0
            cell_diameter = int(self.cell_diameter / self.voxel_length)
            for iteration in range(10):
                radius = 10
                length = 0
                while length <= radius:
                    xi = np.random.randint(0, self.dim.x - 2 * cell_diameter)
                    yi = np.random.randint(0, self.dim.y - 2 * cell_diameter)
                    length = np.sqrt((self.dim.x // 2 - xi) ** 2 + (self.dim.y // 2 - yi) ** 2)
                for x in range(xi, xi + int(cell_diameter)):
                    for y in range(yi, yi + int(cell_diameter)):
                        cell = self.cell_field[x, y, 1]
                        if cell:
                            open_space = False
                            break
                if open_space:
                    concentration_iteration = self.seeding_field[xi, yi, 1]
                    if concentration_iteration >= viral_concentration:
                        viral_concentration = concentration_iteration
                        x_seed = xi
                        y_seed = yi
            if open_space:
                cell = self.new_cell(self.immune_type_id)
                self.cell_field[x_seed:x_seed + int(cell_diameter), y_seed:y_seed + int(cell_diameter), 1] = cell

    def set_seeding_field_name(self, _name: str):
        """
        Set the name of the seeding field

        :param _name: seeding field name
        :return: None
        """
        self._seeding_field_name = _name

    @property
    def seeding_field_name(self):
        """
        Name of the seeding field
        """
        return self._seeding_field_name

    @seeding_field_name.setter
    def seeding_field_name(self, _name: str):
        self.set_seeding_field_name(_name)

    @property
    def seeding_field(self):
        """
        Reference to the seeding field
        """
        return getattr(self.field, self._seeding_field_name)

    def get_immune_seeding_prob(self):
        """
        Get the current probability of seeding an immune cell

        :return: seeding probability
        """
        if self.ir_steppable is not None:
            return self.ir_steppable.get_immune_seeding_prob()
        return 0.0

    def get_immune_removal_prob(self):
        """
        Get the current probability of removing an immune cell

        :return: removal probability
        """
        if self.ir_steppable is not None:
            return self.ir_steppable.get_immune_removal_prob()
        return 0.0

    def on_new_cell(self, _new_cell):
        if _new_cell.type == self.immune_type_id:
            _new_cell.targetVolume = self._cell_volume / self.voxel_length ** 2
            _new_cell.lambdaVolume = self._volume_lm

    def set_cell_volume(self, _val: float):
        """
        Set target volume for newly created cells

        :param _val: target volume for newly created cells
        :return: None
        """
        if _val < 0:
            raise ValueError("Cell target volume must be non-negative")
        self._cell_volume = _val

    @property
    def cell_volume(self):
        """
        Target volume for newly created cells
        """
        return self._cell_volume

    @cell_volume.setter
    def cell_volume(self, _val: float):
        self.set_cell_volume(_val)

    @property
    def cell_diameter(self):
        """
        Approximate cell diameter
        """
        return math.sqrt(self._cell_volume)

    def set_volume_parameter(self, _val: float):
        """
        Set the volume parameter of newly created cells

        :param _val: Volume parameter of newly created cells
        :return: None
        """
        if _val < 0:
            raise ValueError("Volume parameter must be positive")
        self._volume_lm = _val

    @property
    def volume_parameter(self):
        """
        Volume parameter of newly created cells
        """
        return self._volume_lm

    @volume_parameter.setter
    def volume_parameter(self, _val: float):
        self.set_volume_parameter(_val)

    def set_initial_immune_seeding(self, _val: int):
        """
        Set initial immune cell seeding value

        :param _val: Initial immune cell seeding value
        :return: None
        """
        if _val < 0:
            raise ValueError('Initial immune cell seeding must be non-negative')
        self._initial_immune_seeding = _val

    @property
    def initial_immune_seeding(self):
        """
        Initial immune cell seeding value
        """
        return self._initial_immune_seeding

    @initial_immune_seeding.setter
    def initial_immune_seeding(self, _val: int):
        self.set_initial_immune_seeding(_val)


class SimDataSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Plots/writes simulation data of interest

    The frequency of plotting population data in Player can be set with the attribute plot_pop_data_freq. Plotting
    is disabled when plot_pop_data_freq is set to zero.

    The frequency of writing population data to file can be set with the attribute write_pop_data_freq.
    Data is written to the cc3d output directory with the name 'pop_data.dat'.
    Writing is disabled when write_pop_data_freq is set to zero.

    The frequency of plotting diffusion field data in Player can be set with the attribute `plot_med_diff_data_freq`.
    Plotting is disabled when `plot_med_diff_data_freq` is set to zero.

    The frequency of writing diffusion field data to file can be set with the attribute `write_med_diff_data_freq`.
    Data is written to the cc3d output directory with the name 'med_diff_data.dat'.
    Writing is disabled when `write_med_diff_data_freq` is set to zero.

    The frequency of plotting immune response data in Player can be set with the attribute `plot_ir_data_freq`.
    Plotting is disabled when `plot_ir_data_freq` is set to zero.

    The frequency of writing immune response data to file can be set with the attribute `write_ir_data_freq`.
    Data is written to the cc3d output directory with the name 'ir_data.dat'.
    Writing is disabled when `write_ir_data_freq` is set to zero.

    The frequency of plotting death data in Player can be set with the attribute `plot_death_data_freq`.
    Plotting is disabled when `plot_death_data_freq` is set to zero.

    The frequency of writing death data to file can be set with the attribute `write_death_data_freq`.
    Data is written to the cc3d output directory with the name 'death_data.dat'.
    Writing is disabled when `write_death_data_freq` is set to zero.

    The name of each tracked cell type and the virus field can be set with the following attributes
        - uninfected: `uninfected_type_name`
        - infected: `infected_type_name`
        - virus releasing: `virus_releasing_type_name`
        - dead: `dead_type_name`
        - immune: `immune_type_name`
        - virus field: `virus_field_name`
        - cytokine field: `cytokine_field_name`
        - oxidative agent field: `oxidator_field_name`

    Tracking immune response data is performed by calling `get_state_variable_val` on the attribute `ir_steppable`,
    which is a reference to a running recruitment model (e.g., `ImmuneRecruitmentSteppable`).
    `ir_steppable` is automatically retrieved by looking for `ImmuneRecruitmentSteppable.unique_key` in the shared
    steppable dictionary.

    All data writing is performed every quarter-simulation period. All data files are named with the prefix
    module_prefix defined for this module (e.g., 'vivtm_pop_data.dat').
    """

    unique_key = ViralInfectionVTMLib.simdata_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.pop_data_win = None
        self.pop_data_path = None
        self.pop_data = dict()

        self.med_diff_data_win = None
        self.med_diff_data_path = None
        self.med_diff_data = dict()

        self.ir_data_win = None
        self.ir_data_path = None
        self.ir_data = dict()

        self.death_data_win = None
        self.death_data_path = None
        self.death_data = dict()

        self._plot_pop_data_freq = 0
        self._write_pop_data_freq = 0
        self._plot_med_diff_data_freq = 0
        self._write_med_diff_data_freq = 0
        self._plot_ir_data_freq = 0
        self._write_ir_data_freq = 0
        self._plot_death_data_freq = 0
        self._write_death_data_freq = 0

        self.med_diff_key = "MedDiff"
        self.ir_key = "ImmuneResp"
        self.ir_steppable = None

        self.plot_pop_data = False
        self.write_pop_data = False
        self.plot_med_diff_data = False
        self.write_med_diff_data = False
        self.plot_ir_data = False
        self.write_ir_data = False
        self.plot_death_data = False
        self.write_death_data = False

        # Cell death mechanism tracking
        self._death_mech = {'viral': 0,
                            'oxi': 0,
                            'contact': 0,
                            'bystander': 0}

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

        # Initialize defaults
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.set_immune_type_name(immune_type_name)
        self.set_virus_field_name(MainSteppables.virus_field_name)
        self.set_cytokine_field_name(cytokine_field_name)
        self.set_oxidator_field_name(oxidator_field_name)

        self.set_plot_pop_data_freq(ViralInfectionVTMModelInputs.plot_pop_data_freq)
        self.set_write_pop_data_freq(ViralInfectionVTMModelInputs.write_pop_data_freq)
        self.set_plot_med_diff_data_freq(ViralInfectionVTMModelInputs.plot_med_diff_data_freq)
        self.set_write_med_diff_data_freq(ViralInfectionVTMModelInputs.write_med_diff_data_freq)
        self.set_plot_ir_data_freq(ViralInfectionVTMModelInputs.plot_ir_data_freq)
        self.set_write_ir_data_freq(ViralInfectionVTMModelInputs.write_ir_data_freq)
        self.set_plot_death_data_freq(ViralInfectionVTMModelInputs.plot_death_data_freq)
        self.set_write_death_data_freq(ViralInfectionVTMModelInputs.write_death_data_freq)

    def start(self):
        """
        Called once to initialize simulation
        """

        self.plot_pop_data = self._plot_pop_data_freq > 0
        self.write_pop_data = self._write_pop_data_freq > 0
        self.plot_med_diff_data = self._plot_med_diff_data_freq > 0
        self.write_med_diff_data = self._write_med_diff_data_freq > 0
        self.plot_ir_data = self._plot_ir_data_freq > 0
        self.write_ir_data = self._write_ir_data_freq > 0
        self.plot_death_data = self._plot_death_data_freq > 0
        self.write_death_data = self._write_death_data_freq > 0

        # Initialize population data plot if requested
        if self.plot_pop_data:
            self.pop_data_win = self.add_new_plot_window(title='Population data',
                                                         x_axis_title='MCS',
                                                         y_axis_title='Numer of cells',
                                                         x_scale_type='linear',
                                                         y_scale_type='log',
                                                         grid=True,
                                                         config_options={'legend': True})

            self.pop_data_win.add_plot("Uninfected", style='Dots', color='blue', size=5)
            self.pop_data_win.add_plot("Infected", style='Dots', color='red', size=5)
            self.pop_data_win.add_plot("VirusReleasing", style='Dots', color='green', size=5)
            self.pop_data_win.add_plot("Dying", style='Dots', color='yellow', size=5)
            self.pop_data_win.add_plot("ImmuneCell", style='Dots', color='white', size=5)
            self.pop_data_win.add_plot("ImmuneCellActivated", style='Dots', color='purple', size=5)

        if self.plot_med_diff_data:
            self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                              x_axis_title='MCS',
                                                              y_axis_title='Number of diffusive species per volume',
                                                              x_scale_type='linear',
                                                              y_scale_type='log',
                                                              grid=True,
                                                              config_options={'legend': True})

            self.med_diff_data_win.add_plot("MedViral", style='Dots', color='red', size=5)
            self.med_diff_data_win.add_plot("MedCyt", style='Dots', color='blue', size=5)
            self.med_diff_data_win.add_plot("MedOxi", style='Dots', color='green', size=5)

        if self.plot_ir_data:
            self.ir_data_win = self.add_new_plot_window(title='Immune Response Model',
                                                        x_axis_title='MCS',
                                                        y_axis_title='State variable S',
                                                        x_scale_type='linear',
                                                        y_scale_type='linear',
                                                        grid=True)

            self.ir_data_win.add_plot(self.ir_key, style='Dots', color='red', size=5)

        if self.plot_death_data:
            self.death_data_win = self.add_new_plot_window(title='Death data',
                                                           x_axis_title='MCS',
                                                           y_axis_title='Numer of cells',
                                                           x_scale_type='linear',
                                                           y_scale_type='log',
                                                           grid=True,
                                                           config_options={'legend': True})

            self.death_data_win.add_plot("Viral", style='Dots', color='blue', size=5)
            self.death_data_win.add_plot("OxiField", style='Dots', color='red', size=5)
            self.death_data_win.add_plot("Contact", style='Dots', color='green', size=5)
            self.death_data_win.add_plot("Bystander", style='Dots', color='yellow', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath(ViralInfectionVTMLib.module_prefix + 'pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

            if self.write_med_diff_data:
                self.med_diff_data_path = Path(self.output_dir).joinpath(
                    ViralInfectionVTMLib.module_prefix + 'med_diff_data.dat')
                with open(self.med_diff_data_path, 'w'):
                    pass

            if self.write_ir_data:
                self.ir_data_path = Path(self.output_dir).joinpath(ViralInfectionVTMLib.module_prefix + 'ir_data.dat')
                with open(self.ir_data_path, 'w'):
                    pass

            if self.write_death_data:
                self.death_data_path = Path(self.output_dir).joinpath(
                    ViralInfectionVTMLib.module_prefix + 'death_data.dat')
                with open(self.death_data_path, 'w'):
                    pass

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """

        plot_pop_data = self.plot_pop_data and mcs % self._plot_pop_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % self._plot_med_diff_data_freq == 0
        plot_ir_data = self.plot_ir_data and mcs % self._plot_ir_data_freq == 0
        plot_death_data = self.plot_death_data and mcs % self._plot_death_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % self._write_pop_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % self._write_med_diff_data_freq == 0
            write_ir_data = self.write_ir_data and mcs % self._write_ir_data_freq == 0
            write_death_data = self.write_death_data and mcs % self._write_death_data_freq == 0
        else:
            write_pop_data = False
            write_med_diff_data = False
            write_ir_data = False
            write_death_data = False

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.uninfected_type_id))
            num_cells_infected = len(self.cell_list_by_type(self.infected_type_id))
            num_cells_virusreleasing = len(self.cell_list_by_type(self.virus_releasing_type_id))
            num_cells_dying = len(self.cell_list_by_type(self.dead_type_id))
            num_cells_immune = len(self.cell_list_by_type(self.immune_type_id))
            num_cells_immune_act = len([c for c in self.cell_list_by_type(self.immune_type_id)
                                        if c.dict[ViralInfectionVTMLib.activated_cellg_key]])

            # Plot population data plot if requested
            if plot_pop_data:
                if num_cells_uninfected > 0:
                    self.pop_data_win.add_data_point('Uninfected', mcs, num_cells_uninfected)
                if num_cells_infected > 0:
                    self.pop_data_win.add_data_point('Infected', mcs, num_cells_infected)
                if num_cells_virusreleasing > 0:
                    self.pop_data_win.add_data_point('VirusReleasing', mcs, num_cells_virusreleasing)
                if num_cells_dying > 0:
                    self.pop_data_win.add_data_point('Dying', mcs, num_cells_dying)
                if num_cells_immune > 0:
                    self.pop_data_win.add_data_point('ImmuneCell', mcs, num_cells_immune)
                if num_cells_immune_act > 0:
                    self.pop_data_win.add_data_point('ImmuneCellActivated', mcs, num_cells_immune_act)

            # Write population data to file if requested
            if write_pop_data:
                self.pop_data[mcs] = [num_cells_uninfected,
                                      num_cells_infected,
                                      num_cells_virusreleasing,
                                      num_cells_dying,
                                      num_cells_immune,
                                      num_cells_immune_act]

        if plot_med_diff_data or write_med_diff_data:

            # Gather total diffusive amounts
            try:
                med_viral_total = self.get_field_secretor(self._virus_field_name).totalFieldIntegral()
                med_cyt_total = self.get_field_secretor(self._cytokine_field_name).totalFieldIntegral()
                med_oxi_total = self.get_field_secretor(self._oxidator_field_name).totalFieldIntegral()
            except AttributeError:  # Pre-v4.2.1 CC3D
                med_viral_total = 0.0
                med_cyt_total = 0.0
                med_oxi_total = 0.0
                for x, y, z in self.every_pixel():
                    med_viral_total += self.virus_field[x, y, z]
                    med_cyt_total += self.cytokine_field[x, y, z]
                    med_oxi_total += self.oxidator_field[x, y, z]

            # Plot total diffusive viral amount if requested
            if plot_med_diff_data:
                if med_viral_total > 0:
                    self.med_diff_data_win.add_data_point("MedViral", mcs, med_viral_total)
                if med_cyt_total > 0:
                    self.med_diff_data_win.add_data_point("MedCyt", mcs, med_cyt_total)
                if med_oxi_total > 0:
                    self.med_diff_data_win.add_data_point("MedOxi", mcs, med_oxi_total)

            # Write total diffusive viral amount if requested
            if write_med_diff_data:
                self.med_diff_data[mcs] = [med_viral_total,
                                           med_cyt_total,
                                           med_oxi_total]

        if plot_ir_data or write_ir_data:
            if self.ir_steppable is None:
                self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[
                    ImmuneRecruitmentSteppable.unique_key]

            s_val = self.ir_steppable.get_state_variable_val()

            # Plot state variable S if requested
            if plot_ir_data:
                self.ir_data_win.add_data_point(self.ir_key, mcs, s_val)

            # Write state variable S if requested
            if write_ir_data:
                self.ir_data[mcs] = [s_val]

        if plot_death_data or write_death_data:
            num_viral = self._death_mech['viral']
            num_oxi = self._death_mech['oxi']
            num_contact = self._death_mech['contact']
            num_bystander = self._death_mech['bystander']

            # Plot death data if requested
            if plot_death_data:
                if num_viral > 0:
                    self.death_data_win.add_data_point("Viral", mcs, num_viral)
                if num_oxi > 0:
                    self.death_data_win.add_data_point("OxiField", mcs, num_oxi)
                if num_contact > 0:
                    self.death_data_win.add_data_point("Contact", mcs, num_contact)
                if num_bystander > 0:
                    self.death_data_win.add_data_point("Bystander", mcs, num_bystander)

            # Write death data if requested
            if write_death_data:
                self.death_data[mcs] = [num_viral,
                                        num_oxi,
                                        num_contact,
                                        num_bystander]

        # Flush outputs at quarter simulation lengths
        if mcs >= int(self.simulator.getNumSteps() / 4 * self.__flush_counter):
            self.flush_stored_outputs()
            self.__flush_counter += 1

    def on_stop(self):
        """
        Called once when simulation is terminated before completions
        """
        self.finish()

    def finish(self):
        """
        Called once when simulation completes
        """
        self.flush_stored_outputs()

    def data_output_string(self, _data: dict):
        """
        Generate string for data output to file from data dictionary

        :param _data: data dictionary; keys are steps, values are lists of data
        :return: output string to write to file
        """
        mcs_list = list(_data.keys())
        mcs_list.sort()
        f_str = ''
        for mcs in mcs_list:
            f_str += f'{mcs}'
            for v in _data[mcs]:
                f_str += f', {v}'
            f_str += '\n'
        return f_str

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage

        :return: None
        """
        # Each tuple contains the necessary information for writing a set of data to file
        #   1. Boolean for whether we're writing to file at all
        #   2. The path to write the data to
        #   3. The data to write
        output_info = [(self.write_pop_data, self.pop_data_path, self.pop_data),
                       (self.write_med_diff_data, self.med_diff_data_path, self.med_diff_data),
                       (self.write_ir_data, self.ir_data_path, self.ir_data),
                       (self.write_death_data, self.death_data_path, self.death_data)]
        for write_data, data_path, data in output_info:
            if write_data:
                with open(data_path, 'a') as fout:
                    fout.write(self.data_output_string(data))
                    data.clear()

    def track_death_viral(self):
        """
        Increment death count for virally-induced apoptosis
        """
        self._death_mech['viral'] += 1

    def track_death_oxi_field(self):
        """
        Increment death count for oxidative killing
        """
        self._death_mech['oxi'] += 1

    def track_death_contact(self):
        """
        Increment death count for contact killing
        """
        self._death_mech['contact'] += 1

    def track_death_bystander(self):
        """
        Increment death count by the bystander effect
        """
        self._death_mech['bystander'] += 1

    @property
    def death_mech_data(self) -> dict:
        """
        Copy of death mechanism count data
        """
        return self._death_mech.copy()

    def set_plot_pop_data_freq(self, _val):
        """
        Set frequency of plotting population data

        :param _val: Frequency of plotting population data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._plot_pop_data_freq = _val

    @property
    def plot_pop_data_freq(self):
        """
        Frequency of plotting population data
        """
        return self._plot_pop_data_freq

    @plot_pop_data_freq.setter
    def plot_pop_data_freq(self, _val: int):
        self.set_plot_pop_data_freq(_val)

    def set_write_pop_data_freq(self, _val: int):
        """
        Set frequency of writing population data

        :param _val: Frequency of writing population data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._write_pop_data_freq = _val

    @property
    def write_pop_data_freq(self):
        """
        Frequency of writing population data
        """
        return self._write_pop_data_freq

    @write_pop_data_freq.setter
    def write_pop_data_freq(self, _val: int):
        self.set_write_pop_data_freq(_val)

    def set_plot_med_diff_data_freq(self, _val: int):
        """
        Set frequency of plotting diffusion field data

        :param _val: Frequency of plotting diffusion field data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._plot_med_diff_data_freq = _val

    @property
    def plot_med_diff_data_freq(self):
        """
        Frequency of plotting diffusion field data
        """
        return self._plot_med_diff_data_freq

    @plot_med_diff_data_freq.setter
    def plot_med_diff_data_freq(self, _val: int):
        self.set_plot_med_diff_data_freq(_val)

    def set_write_med_diff_data_freq(self, _val):
        """
        Set frequency of writing diffusion field data

        :param _val: Frequency of writing diffusion field data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._write_med_diff_data_freq = _val

    @property
    def write_med_diff_data_freq(self):
        """
        Frequency of writing diffusion field data
        """
        return self._write_med_diff_data_freq

    @write_med_diff_data_freq.setter
    def write_med_diff_data_freq(self, _val: int):
        self.set_write_med_diff_data_freq(_val)

    def set_plot_ir_data_freq(self, _val: int):
        """
        Set frequency of plotting immune response model data

        :param _val: Frequency of plotting immune response model data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._plot_ir_data_freq = _val

    @property
    def plot_ir_data_freq(self):
        """
        Frequency of plotting immune response model data
        """
        return self._plot_ir_data_freq

    @plot_ir_data_freq.setter
    def plot_ir_data_freq(self, _val: int):
        self.set_plot_ir_data_freq(_val)

    def set_write_ir_data_freq(self, _val: int):
        """
        Set frequency of writing immune response model data

        :param _val: Frequency of writing immune response model data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._write_ir_data_freq = _val

    @property
    def write_ir_data_freq(self):
        """
        Frequency of writing immune response model data
        """
        return self._write_ir_data_freq

    @write_ir_data_freq.setter
    def write_ir_data_freq(self, _val: int):
        self.set_write_ir_data_freq(_val)

    def set_plot_death_data_freq(self, _val: int):
        """
        Set frequency of plotting death data

        :param _val: Frequency of plotting death data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._plot_death_data_freq = _val

    @property
    def plot_death_data_freq(self):
        """
        Frequency of plotting death data
        """
        return self._plot_death_data_freq

    @plot_death_data_freq.setter
    def plot_death_data_freq(self, _val: int):
        self.set_plot_death_data_freq(_val)

    def set_write_death_data_freq(self, _val: int):
        """
        Set frequency of writing death data

        :param _val: Frequency of writing death data
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be non-negative")
        self._write_death_data_freq = _val

    @property
    def write_death_data_freq(self):
        """
        Frequency of writing death data
        """
        return self._write_death_data_freq

    @write_death_data_freq.setter
    def write_death_data_freq(self, _val: int):
        self.set_write_death_data_freq(_val)


class CytokineProductionAbsorptionSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements cytokine production/secretion and immune cell activation module, and sets cytokine field data.

    All immune cell types are initialized with the following key: value pairs according to the keys defined in
    `ViralInfectionVTMLib`,
        - `activated_cellg_key`: a flag signifying whether the cell is activated; initialized as False (not activated).
        - `tot_ck_upt_cellg_key`: the current total remembered cytokine uptake of the immune cell; initialized as 0.
        - `ck_production_cellg_key`: the maximum cytokine production value;
          initialized with the settable steppable attribute max_ck_production_immune.
        - `ck_consumption_cellg_key`: the maximum cytokine consumption value;
          initialized with the settable steppable max_ck_consumption_immune.

    Additonally, all infected cells are initialized with the key `ck_production_cellg_key` defined in
    `ViralInfectionVTMLib` and value equal to the settable steppable attribute `max_ck_production_infected`,
    which is the maximum cytokine production of the cell.

    The name of the immune cell type can be set with the attribute `immune_type_name`.

    Requires the module viral replication module as managed by `ViralReplicationSteppable`,
    or one of the same name that defines the variables defined by `vrm_uptake` and `vrm_assembled` in
    `ViralInfectionVTMLib`.

    This steppable uses a callback per registered cell type that is called on each cell of the registered type every
    time step to calculate interactions with the module cytokine field.

    Every secretion callback has the signature `val = secr_func(_steppable, _cell, _mcs)`
        - `_steppable` (CytokineProductionAbsorptionSteppable): self
        - `_cell` (cc3d.cpp.CompuCell.CellG): a cell
        - `_mcs` (int): the current simulation step
        - `val` (float): total produced cytokine by the cell

    Callbacks can be registered for a cell type with `register_secretion_by_type`, and unregistered for a cell type
    with `unregister_secretion_by_type`.

    The step size of the viral replication model per simulation step can be set with the attribute `step_size`.
    This is used as a reference value for calculations, and does not set the step size used in the actual viral
    replication model.

    Total cytokine production can be automatically reported by calling `increment_total_cytokine_count` on the
    attribute `ir_steppable` and passing the increment in total cytokine if `ir_steppable` is set.
    If `ir_steppable` is not set, then it is automatically retrieved by looking for
    `ImmuneRecruitmentSteppable.unique_key` in the shared steppable dictionary.
    If `ir_steppable` is not found in the shared dictionary, then this functionality is ignored.

    By default, requires a diffusion field with module cytokine field name and CC3DML ids 'cytokine_dc' and
    'cytokine_decay' for the field diffusion and decay constants, respectively. These can be set with `set_field_data`.
    Otherwise, these values can be set with the attributes `diffusion_coefficient` and `decay_coefficient`, in which
    case the CC3DML id constraints are no longer relevant.

    By default, module infected and virus-releasing types are assigned the callback `secr_func_infected` and module
    immune cell type is assigned the callback `secr_func_immune`. The immune cell type callback implements the model
    of immune cell activation.

    Model parameters used in secr_func_immune can be set with the attributes as follows
        - `ck_memory_immune`: Immune cell bound cytokine memory parameter
        - `ec50_ck_immune_activation`: Immune cell activation EC50 parameter
        - `minimum_activated_time`: Immune cell minimum activated time
        - `ec50_immune_ck_production`: Immune cell cytokine production EC50 parameter

    Model parameters used in secr_func_infected can be set with the attributes as follows
        - `max_ck_production_infected`: Maximum cytokine production by infected cells

    To track the activated state of immune cells, set `track_model_variables` to True.
    """

    unique_key = ViralInfectionVTMLib.cytokine_secretion_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        self._cytokine_diffusion_id = ''
        self._cytokine_decay_id = ''
        self._secretors_by_type = {}
        self._diffusion_coefficient = None
        self._decay_coefficient = None
        self._step_size = 1.0
        self._ck_memory_immune = 0.0
        self._ec50_ck_immune = 0
        self._ec50_immune_ck_prod = 0
        self._ec50_infecte_ck_prod = 0
        self._minimum_activated_time = 0
        self._max_ck_secrete_infect = 0
        self._max_ck_secrete_im = 0
        self._max_ck_consume = 0
        self._track_model_variables = False

        # Reference to ImmuneResponseSteppable
        self.ir_steppable = None

        # Initialize defaults
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_immune_type_name(immune_type_name)
        self.set_field_data(field_name=cytokine_field_name, diffusion='cytokine_dc', decay='cytokine_decay')
        self.set_step_size(ViralInfectionVTMModelInputs.vr_step_size)
        self.set_ck_memory_immune(ViralInfectionVTMModelInputs.ck_memory_immune)
        self.set_ec50_ck_immune_activation(ViralInfectionVTMModelInputs.EC50_ck_immune)
        self.set_ec50_immune_ck_production(ViralInfectionVTMModelInputs.ec50_immune_ck_prod)
        self.set_ec50_infected_ck_production(ViralInfectionVTMModelInputs.ec50_infecte_ck_prod)
        self.set_minimum_activated_time(ViralInfectionVTMModelInputs.minimum_activated_time)
        self.set_max_ck_production_infected(ViralInfectionVTMModelInputs.max_ck_secrete_infect)
        self.set_max_ck_production_immune(ViralInfectionVTMModelInputs.max_ck_secrete_im)
        self.set_max_ck_consumption_immune(ViralInfectionVTMModelInputs.max_ck_consume)
        self.set_track_model_variables(ViralInfectionVTMModelInputs.track_model_variables)
        self.register_secretion_by_type(MainSteppables.infected_type_name, self.secr_func_infected)
        self.register_secretion_by_type(MainSteppables.virus_releasing_type_name, self.secr_func_infected)
        self.register_secretion_by_type(immune_type_name, self.secr_func_immune)

    def start(self):
        """
        Called once to initialize simulation
        """
        if self._track_model_variables:
            attribute_name = ViralInfectionVTMLib.activated_cellg_key
            field_name = attribute_name.replace(ViralInfectionVTMLib.module_prefix, '')
            self.track_cell_level_scalar_attribute(field_name=field_name, attribute_name=attribute_name)

        # cytokine diff parameters
        diff_factor = self.step_period / (self.voxel_length * self.voxel_length)
        decay_factor = self.step_period
        if self._diffusion_coefficient is None:
            self.get_xml_element(self._cytokine_diffusion_id).cdata = \
                ViralInfectionVTMModelInputs.cytokine_dc * diff_factor
        else:
            self.get_xml_element(self._cytokine_diffusion_id).cdata = self._diffusion_coefficient * diff_factor
        if self._decay_coefficient is None:
            self.get_xml_element(self._cytokine_decay_id).cdata = \
                ViralInfectionVTMModelInputs.cytokine_field_decay * decay_factor
        else:
            self.get_xml_element(self._cytokine_decay_id).cdata = self._decay_coefficient * decay_factor

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = \
                self.shared_steppable_vars[ImmuneRecruitmentSteppable.unique_key]

        # Track the total amount added and subtracted to the cytokine field
        total_ck_inc = 0.0

        for type_name, secr_func in self._secretors_by_type.items():
            for cell in self.cell_list_by_type(getattr(self, type_name.upper())):
                total_ck_inc += secr_func(self, cell, mcs)

        if self.ir_steppable is not None:
            self.ir_steppable.increment_total_cytokine_count(total_ck_inc)

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        """
        Set diffusion field data for cytokine field

        :param field_name: name of the cytokine field (optional)
        :param diffusion: cc3dml id of the diffusion field diffusion coefficient (optional)
        :param decay: cc3dml id of the diffusion field decay coefficient (optional)
        :return: None
        """
        if field_name is not None:
            self.set_cytokine_field_name(field_name)
        if diffusion is not None:
            self._cytokine_diffusion_id = diffusion
        if decay is not None:
            self._cytokine_decay_id = decay

    @staticmethod
    def secr_func_infected(self, _cell, _mcs):
        """
        Secretion function callback for cells registered as infected types by default
        """
        ck_secretor = self.get_field_secretor(self._cytokine_field_name)
        viral_load = ViralInfectionVTMLib.get_assembled_viral_load_inside_cell(_cell, self._step_size)
        produced = _cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] * nCoVUtils.hill_equation(
            viral_load, self.ec50_infected_ck_production, 2)
        res = ck_secretor.secreteInsideCellTotalCount(_cell, produced / _cell.volume)
        return res.tot_amount

    @staticmethod
    def secr_func_immune(self, _cell, _mcs):
        """
        Secretion function callback for cells registered as immune type by default
        """
        total_ck_inc = 0.0
        ck_secretor = self.get_field_secretor(self._cytokine_field_name)

        up_res = ck_secretor.uptakeInsideCellTotalCount(
            _cell, _cell.dict[ViralInfectionVTMLib.ck_consumption_cellg_key] / _cell.volume, 0.1)
        # decay seen ck
        _cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] *= self.ck_memory_immune * self.step_period

        # uptake ck

        # from POV of secretion uptake is negative
        _cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] -= up_res.tot_amount
        total_ck_inc += up_res.tot_amount
        ec50_ck_immune = self.ec50_ck_immune_activation * (
                self.voxel_length ** 3 * ViralInfectionVTMModelInputs.mol_p_L_2_mol_p_um3 *
                ViralInfectionVTMModelInputs.pmol_to_cc3d_au)
        p_activate = nCoVUtils.hill_equation(_cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key],
                                             ec50_ck_immune,
                                             2)

        if rng.uniform() < p_activate and not _cell.dict[ViralInfectionVTMLib.activated_cellg_key]:

            _cell.dict[ViralInfectionVTMLib.activated_cellg_key] = True
            _cell.dict[ViralInfectionVTMLib.time_activation_cellg_key] = _mcs
        elif (_cell.dict[ViralInfectionVTMLib.activated_cellg_key]
              and _mcs - _cell.dict[ViralInfectionVTMLib.time_activation_cellg_key] >
              self._minimum_activated_time / self.step_period):
            _cell.dict[ViralInfectionVTMLib.activated_cellg_key] = False
            _cell.dict[ViralInfectionVTMLib.time_activation_cellg_key] = - 99

        if _cell.dict[ViralInfectionVTMLib.activated_cellg_key]:
            seen_field = self.total_seen_field(self.field.cytokine, _cell)
            ec50_immune_ck_prod = self.ec50_immune_ck_production * (self.voxel_length ** 3 *
                                                                    ViralInfectionVTMModelInputs.mol_p_L_2_mol_p_um3 *
                                                                    ViralInfectionVTMModelInputs.pmol_to_cc3d_au)
            produced = _cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] * nCoVUtils.hill_equation(
                seen_field, ec50_immune_ck_prod, 1)
            sec_res = ck_secretor.secreteInsideCellTotalCount(_cell, produced / _cell.volume)

            total_ck_inc += sec_res.tot_amount

        return total_ck_inc

    def register_secretion_by_type(self, _type_name: str, _secr_func):
        """
        Register a secretion callback by cell type

        :param _type_name: name of cell type
        :param _secr_func: secretion callback
        :return: None
        """
        self._secretors_by_type[_type_name] = _secr_func

    def unregister_secretion_by_type(self, _type_name: str):
        """
        Unregister a secretion callback by cell type

        :param _type_name: name of cell type
        :return: callback of registered cell type
        """
        return self._secretors_by_type.pop(_type_name)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback.

        If a new cell is of registered infected, virus-releasing or immune types, then its cytokine production
        parameter is set.

        If a new cell is of registered immune cell type, total cytokine uptake variable is initialized as zero and
        activated flag is initialized as False as necessary
        """
        if _new_cell.type == self.immune_type_id:
            if ViralInfectionVTMLib.tot_ck_upt_cellg_key not in _new_cell.dict.keys():
                _new_cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] = 0
            if ViralInfectionVTMLib.activated_cellg_key not in _new_cell.dict.keys():
                _new_cell.dict[ViralInfectionVTMLib.activated_cellg_key] = False
        self.on_set_cell_type(cell=_new_cell, old_type=None)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback. If a cell's type changes to the registered virus-releasing type, then
        the cytokine production model parameter is updated
        """
        if cell.type in [self.infected_type_id, self.virus_releasing_type_id]:
            # update cytokine param
            cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = self._max_ck_secrete_infect * (
                    self.step_period * ViralInfectionVTMModelInputs.mol_p_L_2_mol_p_um3 *
                    self.voxel_length ** 3 * ViralInfectionVTMModelInputs.pmol_to_cc3d_au)
        elif cell.type == self.immune_type_id:
            conv_fact = (self.step_period * ViralInfectionVTMModelInputs.mol_p_L_2_mol_p_um3 *
                         self.voxel_length ** 3 * ViralInfectionVTMModelInputs.pmol_to_cc3d_au)
            cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = self._max_ck_secrete_im * conv_fact
            cell.dict[ViralInfectionVTMLib.ck_consumption_cellg_key] = self._max_ck_consume * conv_fact

    def set_diffusion_coefficient(self, _val: float):
        """
        Set diffusion coefficient, in units microns^2/s

        :param _val: Diffusion coefficient
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Diffusion coefficient must be positive")
        self._diffusion_coefficient = _val

    @property
    def diffusion_coefficient(self):
        """
        Diffusion coefficient, in units microns^2/s
        """
        return self._diffusion_coefficient

    @diffusion_coefficient.setter
    def diffusion_coefficient(self, _val: float):
        self.set_diffusion_coefficient(_val)

    def set_decay_coefficient(self, _val: float):
        """
        Set decay coefficient, in units 1/s

        :param _val: Decay coefficient
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Decay coefficient must be positive")
        self._decay_coefficient = _val

    @property
    def decay_coefficient(self):
        """
        Decay coefficient, in units 1/s
        """
        return self._decay_coefficient

    @decay_coefficient.setter
    def decay_coefficient(self, _val: float):
        self.set_decay_coefficient(_val)

    def set_step_size(self, _val: float):
        """
        Set viral replication model step size

        :param _val: step size
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Step size must be positive")
        self._step_size = _val

    @property
    def step_size(self):
        """
        Viral replication model step size
        """
        return self._step_size

    @step_size.setter
    def step_size(self, _val: float):
        self.set_step_size(_val)

    def set_ck_memory_immune(self, _val: float):
        """
        Set immune cell bound cytokine memory parameter

        :param _val: Immune cell bound cytokine memory parameter
        :return: None
        """
        if not 0.0 <= _val < 1.0:
            raise ValueError("Value must be in [0, 1)")
        self._ck_memory_immune = _val

    @property
    def ck_memory_immune(self):
        """
        Immune cell bound cytokine memory parameter
        """
        return self._ck_memory_immune

    @ck_memory_immune.setter
    def ck_memory_immune(self, _val: float):
        self.set_ck_memory_immune(_val)

    def set_ec50_ck_immune_activation(self, _val: float):
        """
        Set immune cell activation EC50 parameter, in units pM

        :param _val: Immune cell activation EC50 parameter
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Value must be positive")
        self._ec50_ck_immune = _val

    @property
    def ec50_ck_immune_activation(self):
        """
        Immune cell activation EC50 parameter, in units pM
        """
        return self._ec50_ck_immune

    @ec50_ck_immune_activation.setter
    def ec50_ck_immune_activation(self, _val: float):
        self.set_ec50_ck_immune_activation(_val)

    def set_ec50_immune_ck_production(self, _val: float):
        """
        Set immune cell cytokine production EC50 parameter, in units pM

        :param _val: Immune cell cytokine production EC50 parameter
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._ec50_immune_ck_prod = _val

    @property
    def ec50_immune_ck_production(self):
        """
        Immune cell cytokine production EC50 parameter, in units pM
        """
        return self._ec50_immune_ck_prod

    @ec50_immune_ck_production.setter
    def ec50_immune_ck_production(self, _val: float):
        self.set_ec50_immune_ck_production(_val)

    def set_ec50_infected_ck_production(self, _val: float):
        """
        Set infected cell cytokine production EC50 parameter

        :param _val: Infected cell cytokine production EC50 parameter
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._ec50_infecte_ck_prod = _val

    @property
    def ec50_infected_ck_production(self):
        """
        Infected cell cytokine production EC50 parameter
        """
        return self._ec50_infecte_ck_prod

    @ec50_infected_ck_production.setter
    def ec50_infected_ck_production(self, _val: float):
        self.set_ec50_infected_ck_production(_val)

    def set_minimum_activated_time(self, _val: float):
        """
        Set immune cell minimum activated time, in units s

        :param _val: Immune cell minimum activated time
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be positive")
        self._minimum_activated_time = _val

    @property
    def minimum_activated_time(self):
        """
        Immune cell minimum activated time, in units s
        """
        return self._minimum_activated_time

    @minimum_activated_time.setter
    def minimum_activated_time(self, _val: float):
        self.set_minimum_activated_time(_val)

    def set_max_ck_production_infected(self, _val: float):
        """
        Set maximum cytokine production by infected cells, in units pM/s

        :param _val: Maximum cytokine production by infected cells
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be positive")
        self._max_ck_secrete_infect = _val

    @property
    def max_ck_production_infected(self):
        """
        Maximum cytokine production by infected cells, in units pM/s
        """
        return self._max_ck_secrete_infect

    @max_ck_production_infected.setter
    def max_ck_production_infected(self, _val: float):
        self.set_max_ck_production_infected(_val)

    def set_max_ck_production_immune(self, _val: float):
        """
        Set maximum cytokine production by immune cells, in units pM/s

        :param _val: Maximum cytokine production by immune cells
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be positive")
        self._max_ck_secrete_im = _val

    @property
    def max_ck_production_immune(self):
        """
        Maximum cytokine production by immune cells, in units pM/s
        """
        return self._max_ck_secrete_im

    @max_ck_production_immune.setter
    def max_ck_production_immune(self, _val: float):
        self.set_max_ck_production_immune(_val)

    def set_max_ck_consumption_immune(self, _val: float):
        """
        Set maximum cytokine consumption by immune cells, in units pM/s

        :param _val: Maximum cytokine consumption by immune cells
        :return: None
        """
        if _val < 0:
            raise ValueError("Value must be positive")
        self._max_ck_consume = _val

    @property
    def max_ck_consumption_immune(self):
        """
        Maximum cytokine consumption by immune cells, in units pM/s
        """
        return self._max_ck_consume

    @max_ck_consumption_immune.setter
    def max_ck_consumption_immune(self, _val: float):
        self.set_max_ck_consumption_immune(_val)

    def set_track_model_variables(self, _val: bool):
        """
        Set flag to track model variables in cells

        :param _val: flag value: enables when True
        :return: None
        """
        self._track_model_variables = _val

    @property
    def track_model_variables(self):
        """
        Flag to track model variables in cells
        """
        return self._track_model_variables

    @track_model_variables.setter
    def track_model_variables(self, _val: bool):
        self.set_track_model_variables(_val)


class ImmuneRecruitmentSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell recruitment module

    Note that total cytokine is currently tracked elsewhere by counting uptake and secretion, and by applying the
    field decay rate applied to the cytokine field. This is only relevant for periodic and zero-flux boundary conditions
    for the diffusive cytokine field.

    The step size of the immune recruitment model per simulation step can be set with the attribute `step_size`.

    The probability scaling factor can be set with the attribute `prob_scaling_factor`.

    By default, looks for a CC3DML tag 'cytokine_decay' to get the field decay coefficient of cytokine.
    The decay coefficient can instead be set with the attribute `cytokine_decay`.

    To manage the total cytokine manually, use `increment_total_cytokine_count`.

    The immune recruitment model can be initialized with different model parameters using the following attributes,
        - `add_coeff`: Immune recruitment model add coefficient
        - `subtract_coeff`: Immune recruitment model subtract coefficient
        - `delay_coeff`: Delay coefficient
        - `decay_coeff`: Decay coefficient
        - `transmission_coeff`: Transmission coefficient

    The recruitment model can be changed by setting/overriding `recruitment_model_generator`.

    This steppable interacts with the recruitment model through the following variables,
        - `numImmuneCells`: the current number of immune cells is passed to the recruitment model.
        - `totalCytokine`: the current total cytokine is passed to the recruitment model.
        - `S`: a state variable describing the state of the immune response is retrieved from the recruitment model.

    """

    unique_key = ViralInfectionVTMLib.ir_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to solver
        self.__rr = None

        # Running value of total cytokine; to be updated externally through accessor
        self.__total_cytokine = 0.0

        self.__ck_decay = None
        self._step_size = 1.0
        self._add_coeff = 0.0
        self._subtract_coeff = 0.0
        self._delay_coeff = 0.0
        self._decay_coeff = 0.0
        self._transmission_coeff = 0.0
        self._prob_scaling_factor = 0.0

        # Initialize defaults
        self.set_immune_type_name(immune_type_name)
        self.set_step_size(ViralInfectionVTMModelInputs.vr_step_size)
        self.set_add_coeff(ViralInfectionVTMModelInputs.ir_add_coeff)
        self.set_subtract_coeff(ViralInfectionVTMModelInputs.ir_subtract_coeff)
        self.set_delay_coeff(ViralInfectionVTMModelInputs.ir_delay_coeff)
        self.set_decay_coeff(ViralInfectionVTMModelInputs.ir_decay_coeff)
        self.set_transmission_coeff(ViralInfectionVTMModelInputs.ir_transmission_coeff)
        self.set_prob_scaling_factor(ViralInfectionVTMModelInputs.ir_prob_scaling_factor)

    def start(self):
        """
        Called once to initialize simulation
        """
        if self.__ck_decay is None:
            self.__ck_decay = float(self.get_xml_element('cytokine_decay').cdata)

        # Initialize model
        self.__init_fresh_recruitment_model()

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """

        # Update total count of immune cells
        num_immune_cells = len(self.cell_list_by_type(self.immune_type_id))

        # Apply consumption / transmission decay to running total
        total_cytokine_decayed = self.__total_cytokine * self.__ck_decay
        self.__total_cytokine -= total_cytokine_decayed

        # Update model
        total_cytokine_transmitted = self._transmission_coeff * total_cytokine_decayed
        self.update_running_recruitment_model(num_immune_cells, total_cytokine_transmitted)

    def __init_fresh_recruitment_model(self):
        # Generate solver instance
        self.register_ode_model(model_name=ViralInfectionVTMLib.ir_model_name,
                                model_fcn=self.recruitment_model_generator(),
                                step_size=self._step_size)

        # Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == ViralInfectionVTMLib.ir_model_name:
                self.__rr = rr

    def recruitment_model_generator(self):
        """
        Get the model string generator

        :return: Model string generator
        """
        def model_fcn():
            """
            Callback to the model string generator
            """
            return ViralInfectionVTMLib.immune_recruitment_model_string(
                self._add_coeff * self.step_period,
                self._subtract_coeff * self.step_period,
                self._delay_coeff / self.step_period,
                self._decay_coeff * self.step_period)
        return model_fcn

    def update_running_recruitment_model(self, num_immune_cells, total_cytokine):
        """
        Update running immune recruitment model by one step

        :param num_immune_cells: Number of immune cells
        :param total_cytokine: Total cytokine
        :return: None
        """
        self.__rr['numImmuneCells'] = num_immune_cells
        self.__rr['totalCytokine'] = total_cytokine
        self.timestep_ode_model(model_name=ViralInfectionVTMLib.ir_model_name)

    def get_state_variable_val(self):
        """
        Get immune recruitment model state variable

        :return: Immune recruitment model state variable
        """
        return self.__rr['S']

    def get_immune_seeding_prob(self):
        """
        Returns probability of immune cell seeding due to local and global recruitment

        Probability is only non-zero if the state variable `S` is positive, in which case
        the probability is an error function of `S`

        :return: probability of immune cell seeding due to local and global recruitment
        """
        s_val = self.get_state_variable_val()
        if s_val < 0:
            return 0.0
        else:
            return math.erf(self._prob_scaling_factor * s_val)

    def get_immune_removal_prob(self):
        """
        Returns probability of immune cell removal due to local and global recruitment

        Probability is only non-zero if the state variable `S` is negative, in which case
        the probability is an error function of `S`

        :return: probability of immune cell removal due to local and global recruitment
        """
        s_val = self.get_state_variable_val()
        if s_val > 0:
            return 0.0
        else:
            return math.erf(- self._prob_scaling_factor * s_val)

    def increment_total_cytokine_count(self, _inc_amount):
        """
        Increment total count of cytokine

        :param _inc_amount: total amount increment
        :return: None
        """
        self.__total_cytokine = max(0.0, self.__total_cytokine + _inc_amount)

    def set_add_coeff(self, _val: float):
        """
        Set immune recruitment model add coefficient

        :param _val: Immune recruitment model add coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._add_coeff = _val

    @property
    def add_coeff(self):
        """
        Immune recruitment model add coefficient
        """
        return self._add_coeff

    @add_coeff.setter
    def add_coeff(self, _val: float):
        self.set_add_coeff(_val)

    def set_subtract_coeff(self, _val: float):
        """
        Set immune recruitment model subtract coefficient

        :param _val: Immune recruitment model subtract coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._subtract_coeff = _val

    @property
    def subtract_coeff(self):
        """
        Immune recruitment model subtract coefficient
        """
        return self._subtract_coeff

    @subtract_coeff.setter
    def subtract_coeff(self, _val: float):
        self.set_subtract_coeff(_val)

    def set_delay_coeff(self, _val: float):
        """
        Set immune recruitment model delay coefficient

        :param _val: Immune recruitment model delay coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._delay_coeff = _val

    @property
    def delay_coeff(self):
        """
        Immune recruitment model delay coefficient
        """
        return self._delay_coeff

    @delay_coeff.setter
    def delay_coeff(self, _val: float):
        self.set_delay_coeff(_val)

    def set_decay_coeff(self, _val: float):
        """
        Set immune recruitment model decay coefficient

        :param _val: Immune recruitment model decay coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._decay_coeff = _val

    @property
    def decay_coeff(self):
        """
        Immune recruitment model decay coefficient
        """
        return self._decay_coeff

    @decay_coeff.setter
    def decay_coeff(self, _val: float):
        self.set_decay_coeff(_val)

    def set_transmission_coeff(self, _val: float):
        """
        Set immune recruitment model transmission coefficient

        :param _val: Immune recruitment model transmission coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._transmission_coeff = _val

    @property
    def transmission_coeff(self):
        """
        Immune recruitment model transmission coefficient
        """
        return self._transmission_coeff

    @transmission_coeff.setter
    def transmission_coeff(self, _val: float):
        self.set_transmission_coeff(_val)

    def set_prob_scaling_factor(self, _val: float):
        """
        Set immune recruitment model probability scaling coefficient

        :param _val: Immune recruitment model probability scaling coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Value must be non-negative")
        self._prob_scaling_factor = _val

    @property
    def prob_scaling_factor(self):
        """
        Immune recruitment model probability scaling coefficient
        """
        return self._prob_scaling_factor

    @prob_scaling_factor.setter
    def prob_scaling_factor(self, _val: float):
        self.set_prob_scaling_factor(_val)

    def set_cytokine_decay(self, _val: float):
        """
        Set cytokine decay rate

        :param _val: cytokine decay rate
        :return: None
        """
        if _val < 0:
            raise ValueError('Cytokine decay rate must be non-negative')
        self.__ck_decay = _val

    @property
    def cytokine_decay(self):
        """
        Cytokine decay rate
        """
        return self.__ck_decay

    @cytokine_decay.setter
    def cytokine_decay(self, _val: float):
        self.set_cytokine_decay(_val)

    def set_step_size(self, _val: float):
        """
        Set recruitment model step size

        :param _val: step size
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Step size must be positive")
        self._step_size = _val

    @property
    def step_size(self):
        """
        Immune recruitment model step size
        """
        return self._step_size

    @step_size.setter
    def step_size(self, _val: float):
        self.set_step_size(_val)


class oxidationAgentModelSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell oxidizing agent cytotoxicity module

    This steppable uses a callback per registered cell type that is called on each cell of the registered type every
    time step to calculate interactions with the module oxidator field.

    Every secretion callback has the signature secr_func(_steppable, _cell, _mcs)
    - _steppable (oxidationAgentModelSteppable): self
    - _cell (cc3d.cpp.CompuCell.CellG): a cell
    - _mcs (int): the current simulation step

    Callbacks can be registered for a cell type with register_secretion_by_type, and unregistered for a cell type with
    unregister_secretion_by_type.

    Requires a cytokine diffusion field.
    The name of the cytokine field can be set with the attribute cytokine_field_name.

    By default, requires a diffusion field with module oxidator field name and CC3DML ids 'oxi_dc' and
    'oxi_decay' for the field diffusion and decay constants, respectively. These can be set with set_field_data.
    Otherwise, these values can be set with the attributes diffusion_coefficient and decay_coefficient, in which case
    the CC3DML id constraints are no longer relevant.
    The name of the oxidator field can be set with the attribute oxidator_field_name.

    By default, module uninfected, infected and virus-releasing types are assigned the callback secr_func_epithelial
    and module immune cell type is assigned the callback secr_func_immune.

    secr_func_immune requires a value with key of value activated_cellg_key defined in ViralInfectionVTMLib in the
    dictionary of all cells of the immune type, which is a flag describing the activated state of the immune cell.
    In ordinary usage, this dictionary entry is managed by CytokineProductionAbsorptionSteppable.

    secr_func_immune uses the settable attribute oxi_secr_thr, which is the cytokine threshold above which
    activated immune cells release into the oxidator field.

    secr_func_epithelial uses the settable attribute death_threshold, which is the threshold of oxidative agent above
    which cells of a registered type dies.

    Reports the occurrence of death if the attribute simdata_steppable is set or a module SimDataSteppable instance is
    found in the shared dictionary.
    If simdata_steppable is set, then the referred instance must have the method track_death_oxi_field, which it will
    call without arguments when reporting death.

    If the attribute track_model_variables is set to True, then death by oxidative killing is tracked in each cell
    of a corresponding registered type.
    """

    unique_key = ViralInfectionVTMLib.oxidation_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self._oxidator_diffusion_id = ''
        self._oxidator_decay_id = ''
        self._secretors_by_type = {}
        self._oxi_death_thr = 0.0
        self._secr_rate = 0.0
        self._oxi_secr_thr = 0.0
        self._diffusion_coefficient = None
        self._decay_coefficient = None
        self._track_model_variables = False

        # Initialize default data
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.set_cytokine_field_name(cytokine_field_name)
        self.set_field_data(field_name=oxidator_field_name, diffusion='oxi_dc', decay='oxi_decay')
        self.set_death_threshold(ViralInfectionVTMModelInputs.oxi_death_thr)
        self.set_secr_rate(ViralInfectionVTMModelInputs.max_oxi_secrete)
        self.set_oxi_secr_thr(ViralInfectionVTMModelInputs.oxi_secr_thr)
        self.set_track_model_variables(ViralInfectionVTMModelInputs.track_model_variables)
        self.register_secretion_by_type(MainSteppables.uninfected_type_name, self.secr_func_epithelial)
        self.register_secretion_by_type(MainSteppables.infected_type_name, self.secr_func_epithelial)
        self.register_secretion_by_type(MainSteppables.virus_releasing_type_name, self.secr_func_epithelial)
        self.register_secretion_by_type(immune_type_name, self.secr_func_immune)

    def start(self):
        """
        Called once to initialize simulation
        """
        if self._track_model_variables:
            attribute_name = ViralInfectionVTMLib.oxi_killed_cellg_key
            field_name = attribute_name.replace(ViralInfectionVTMLib.module_prefix, '')
            self.track_cell_level_scalar_attribute(field_name=field_name, attribute_name=attribute_name)

        diff_factor = self.step_period / (self.voxel_length * self.voxel_length)
        decay_factor = self.step_period
        if self._diffusion_coefficient is None:
            self.get_xml_element(self._oxidator_diffusion_id).cdata = ViralInfectionVTMModelInputs.oxi_dc * diff_factor
        else:
            self.get_xml_element(self._oxidator_diffusion_id).cdata = self._diffusion_coefficient * diff_factor
        if self._decay_coefficient is None:
            self.get_xml_element(self._oxidator_decay_id).cdata = ViralInfectionVTMModelInputs.oxi_decay * decay_factor
        else:
            self.get_xml_element(self._oxidator_decay_id).cdata = self._decay_coefficient * decay_factor

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        if self.simdata_steppable is None:
            try:
                self.simdata_steppable: SimDataSteppable = self.shared_steppable_vars[SimDataSteppable.unique_key]
            except KeyError:
                pass

        for type_name, secr_func in self._secretors_by_type.items():
            for cell in self.cell_list_by_type(getattr(self, type_name.upper())):
                secr_func(self, cell, mcs)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback. Manage tracked data.
        """
        if old_type == self.dead_type_id:
            cell.dict[ViralInfectionVTMLib.oxi_killed_cellg_key] = False

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        """
        Set diffusion field data for oxidative field

        :param field_name: name of the oxidative field (optional)
        :param diffusion: cc3dml id of the diffusion field diffusion coefficient (optional)
        :param decay: cc3dml id of the diffusion field decay coefficient (optional)
        :return: None
        """
        if field_name is not None:
            self.set_oxidator_field_name(field_name)
        if diffusion is not None:
            self._oxidator_diffusion_id = diffusion
        if decay is not None:
            self._oxidator_decay_id = decay

    def register_secretion_by_type(self, _type_name: str, _secr_func):
        """
        Register a secretion callback by cell type

        :param _type_name: name of cell type
        :param _secr_func: secretion callback
        :return: None
        """
        self._secretors_by_type[_type_name] = _secr_func

    def unregister_secretion_by_type(self, _type_name: str):
        """
        Unregister a secretion callback by cell type

        :param _type_name: name of cell type
        :return: callback of registered cell type
        """
        return self._secretors_by_type.pop(_type_name)

    @staticmethod
    def secr_func_epithelial(self, _cell, _mcs):
        """
        Secretion function callback for cells registered as epithelial types by default
        """
        seen_field = self.total_seen_field(self.oxidator_field, _cell)
        oxi_death_thr = self._oxi_death_thr * self.voxel_length ** 3 * (
                ViralInfectionVTMModelInputs.mol_p_L_2_mol_p_um3 * ViralInfectionVTMModelInputs.pmol_to_cc3d_au)
        if seen_field >= oxi_death_thr:
            self.kill_cell(cell=_cell)
            _cell.dict[ViralInfectionVTMLib.oxi_killed_cellg_key] = True
            if self.simdata_steppable is not None:
                self.simdata_steppable.track_death_oxi_field()

    @staticmethod
    def secr_func_immune(self, _cell, _mcs):
        """
        Secretion function callback for cells registered as immune type by default
        """
        oxi_secretor = self.get_field_secretor(self._oxidator_field_name)
        if _cell.dict[ViralInfectionVTMLib.activated_cellg_key]:
            seen_field = self.total_seen_field(self.cytokine_field, _cell)
            m_to_au = self.voxel_length ** 3 * (ViralInfectionVTMModelInputs.mol_p_L_2_mol_p_um3 *
                                                ViralInfectionVTMModelInputs.pmol_to_cc3d_au)
            if seen_field > self._oxi_secr_thr * m_to_au:
                secr_rate = self._secr_rate * self.step_period * m_to_au
                oxi_secretor.secreteInsideCellTotalCount(_cell, secr_rate / _cell.volume)

    def set_death_threshold(self, _val: float):
        """
        Set threshold for death by oxidative killing, in units pM

        :param _val: Threshold for death by oxidative killing
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Oxidative agent death threshold must be non-negative")
        self._oxi_death_thr = _val

    @property
    def death_threshold(self):
        """
        Threshold for death by oxidative killing, in units pM
        """
        return self._oxi_death_thr

    @death_threshold.setter
    def death_threshold(self, _val: float):
        self.set_death_threshold(_val)

    def set_secr_rate(self, _val: float):
        """
        Set oxidative field release rate, in units pM/s

        :param _val: Oxidative field release rate
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Oxidative agent secretion rate must be non-negative")
        self._secr_rate = _val

    @property
    def secr_rate(self):
        """
        Oxidative field release rate, in units pM/s
        """
        return self._secr_rate

    @secr_rate.setter
    def secr_rate(self, _val: float):
        self.set_secr_rate(_val)

    def set_oxi_secr_thr(self, _val: float):
        """
        Set cytokine threshold for oxidative release, in units pM

        :param _val: Cytokine threshold for oxidative release
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Oxidative agent threshold must be non-negative")
        self._oxi_secr_thr = _val

    @property
    def oxi_secr_thr(self):
        """
        Cytokine threshold for oxidative release, in units pM
        """
        return self._oxi_secr_thr

    @oxi_secr_thr.setter
    def oxi_secr_thr(self, _val: float):
        self.set_oxi_secr_thr(_val)

    def set_diffusion_coefficient(self, _val: float):
        """
        Set oxidative agent diffusion coefficient

        :param _val: Oxidative agent diffusion coefficient
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Diffusion coefficient must be positive")
        self._diffusion_coefficient = _val

    @property
    def diffusion_coefficient(self):
        """
        Oxidative agent diffusion coefficient
        """
        return self._diffusion_coefficient

    @diffusion_coefficient.setter
    def diffusion_coefficient(self, _val: float):
        self.set_diffusion_coefficient(_val)

    def set_decay_coefficient(self, _val: float):
        """
        Set oxidative agent decay coefficient

        :param _val: Oxidative agent decay coefficient
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Decay coefficient must be positive")
        self._decay_coefficient = _val

    @property
    def decay_coefficient(self):
        """
        Oxidative agent decay coefficient
        """
        return self._decay_coefficient

    @decay_coefficient.setter
    def decay_coefficient(self, _val: float):
        self.set_decay_coefficient(_val)

    def set_track_model_variables(self, _val: bool):
        """
        Set flag to track model variables in cells

        :param _val: flag value: enables when True
        :return: None
        """
        self._track_model_variables = _val

    @property
    def track_model_variables(self):
        """
        Flag to track model variables in cells
        """
        return self._track_model_variables

    @track_model_variables.setter
    def track_model_variables(self, _val: bool):
        self.set_track_model_variables(_val)
