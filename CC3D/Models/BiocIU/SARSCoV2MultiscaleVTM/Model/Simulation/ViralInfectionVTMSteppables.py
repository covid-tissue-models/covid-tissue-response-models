"""
Defines module steppables

Steppables
==========

CellInitializerSteppable
------------------------
Description: Initializes an epithelial sheet with an infection

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import CellInitializerSteppable

steppable = CellInitializerSteppable(frequency=1)

# Change initial infection conditions from default

steppable.random_infected_fraction(0.05)

CompuCellSetup.register_steppable(steppable=steppable)

VirusFieldInitializerSteppable
------------------------------
Description: Initializes virus field data and properties

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import VirusFieldInitializerSteppable

CompuCellSetup.register_steppable(steppable=VirusFieldInitializerSteppable(frequency=1))

ViralInternalizationSteppable
-----------------------------
Description: Performs viral internalization

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import ViralInternalizationSteppable

CompuCellSetup.register_steppable(steppable=ViralInternalizationSteppable(frequency=1))

EclipsePhaseSteppable
---------------------
Description: Performs viral eclipse phase

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import EclipsePhaseSteppable

CompuCellSetup.register_steppable(steppable=EclipsePhaseSteppable(frequency=1))

ViralDeathSteppable
-------------------
Description: Performs viral death

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import ViralDeathSteppable

CompuCellSetup.register_steppable(steppable=ViralDeathSteppable(frequency=1))

ViralReleaseSteppable
---------------------
Description: Performs viral release

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import ViralReleaseSteppable

CompuCellSetup.register_steppable(steppable=ViralReleaseSteppable(frequency=1))

SimDataSteppable
----------------
Description: Plots/writes simulation data of interest

Usage: In ViralInfectionVTM.py, add the following

from ViralInfectionVTMSteppables import SimDataSteppable

CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))
"""

import random

# Import project libraries and classes
import ViralInfectionVTMModelInputs

# Import toolkit
from nCoVToolkit import nCoVSteppableBase, nCoVUtils

# Unique steppable keys
cell_initializer_key = 'std_cell_initializer_steppable'  # CellInitializerSteppable
virus_field_initializer_key = 'std_virus_field_initializer_steppable'  # VirusFieldInitializerSteppable
viral_internalization_key = 'std_viral_internalization_steppable'  # ViralInternalizationSteppable
eclipse_phase_key = 'std_eclipse_phase_steppable'  # EclipsePhaseSteppable
viral_death_key = 'std_viral_death_steppable'  # ViralDeathSteppable
viral_release_key = 'std_viral_release_steppable'  # ViralReleaseSteppable
sim_data_key = 'std_sim_data_steppable'  # SimDataSteppable

# Standard epithelial cell type names
uninfected_type_name = 'Uninfected'
infected_type_name = 'Infected'
virus_releasing_type_name = 'VirusReleasing'
dead_type_name = 'Dying'
# Standard virus field name
virus_field_name = 'Virus'


class CellInitializerSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Initializes an epithelial sheet with an infection

    Currently available initial infection modes are
       - no_initial_infection: no initial infection
       - single_infected_cell: infect a single cell at the center of a sheet
       - random_infected_fraction: infect a randomly selected fraction of an epithelial population
       - random_infected_probability: initially infect each cell with a probability
    """

    unique_key = cell_initializer_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self._initialization_mode = None
        self._infection_mode = None
        self._uninfected_type_name = ''
        self._infected_type_name = ''
        self._cell_diameter = 0
        self._volume_lm = 0.0

        self._aux_infect_vars = None

        # Set default data
        self.set_uninfected_type_name(uninfected_type_name)
        self.set_infected_type_name(infected_type_name)
        self.set_cell_diameter(ViralInfectionVTMModelInputs.cell_diameter)
        self.set_volume_parameter(ViralInfectionVTMModelInputs.volume_lm)
        self.initialize_sheet()
        self.random_infected_fraction(0.01)

    def start(self):
        """
        Called once to initialize simulation
        """
        # Do cell initialization
        self.initialize_cells()

        # Do initial infection
        self.initialize_infection()

    def initialize_cells(self):
        """
        Initialize cells according to selected initialization mode

        :return: None
        """
        if self._initialization_mode is not None:
            self._initialization_mode()

    def initialize_infection(self):
        """
        Initialize infection according to selected initialization mode

        :return: None
        """
        if self._infection_mode is not None:
            self._infection_mode()

    def _initialize_sheet(self):
        # Enforce compatible lattice dimensions with epithelial cell size
        cdiam = int(self.cell_diameter / self.voxel_length)
        assert self.dim.x % cdiam == 0 and self.dim.y % cdiam == 0, \
            f'Lattice dimensions must be multiples of the unitless cell diameter (currently cell_diameter = {cdiam})'

        for x in range(0, self.dim.x, cdiam):
            for y in range(0, self.dim.y, cdiam):
                cell = self.new_cell(self.uninfected_type_id)
                self.cellField[x:x + cdiam, y:y + cdiam, 0] = cell

                cell.targetVolume = self.cell_volume
                cell.lambdaVolume = self.volume_parameter

    def initialize_sheet(self):
        """
        Initialize an epithelial sheet

        :return: None
        """
        self._initialization_mode = self._initialize_sheet

    def _no_initial_infection(self):
        pass

    def _single_infected_cell(self):
        cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
        self.set_cell_type(cell, self.infected_type_id)

    def _random_infected_fraction(self):
        if self._aux_infect_vars is None:
            raise AttributeError("Initial infection fraction not set")
        elif not (0 <= self._aux_infect_vars <= 1):
            raise ValueError("Invalid initial infection fraction: must be in [0, 1]")
        cell_list = [c for c in self.cell_list_by_type(self.uninfected_type_id)]
        num_to_infect = int(len(cell_list) * self._aux_infect_vars)
        random.shuffle(cell_list)
        cells_to_infect = cell_list[0:num_to_infect]
        for cell in cells_to_infect:
            self.set_cell_type(cell, self.infected_type_id)

    def _random_infected_probability(self):
        if self._aux_infect_vars is None:
            raise AttributeError("Initial infection probability not set")
        elif not (0 <= self._aux_infect_vars <= 1):
            raise ValueError("Invalid initial infection probability: must be in [0, 1]")
        for cell in self.cell_list_by_type(self.uninfected_type_id):
            if random.random() <= self._aux_infect_vars:
                self.set_cell_type(cell, self.infected_type_id)

    def no_initial_infection(self):
        """
        Use no initial infection

        :return: None
        """
        self._infection_mode = self._no_initial_infection
        self._aux_infect_vars = None

    def single_infected_cell(self):
        """
        Infect a single cell at the center of a sheet

        :return: None
        """
        self._infection_mode = self._single_infected_cell
        self._aux_infect_vars = None

    def random_infected_fraction(self, _frac: float):
        """
        Infect a randomly selected fraction of an epithelial population

        :param _frac: Fraction of initially infected cells
        :return: None
        """
        if not (0 <= _frac <= 1):
            raise ValueError("Invalid initial infection fraction: must be in [0, 1]")
        self._infection_mode = self._random_infected_fraction
        self._aux_infect_vars = _frac

    def random_infected_probability(self, _prob: float):
        """
        Initially infect each cell with a probability

        :param _prob: Probability of being initially infected
        :return: None
        """
        if not (0 <= _prob <= 1):
            raise ValueError("Invalid initial infection probability: must be in [0, 1]")
        self._infection_mode = self._random_infected_probability
        self._aux_infect_vars = _prob

    @property
    def uninfected_type_id(self) -> int:
        """
        Id of the uninfected cell type according to a cc3d simulation
        """
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def infected_type_id(self) -> int:
        """
        Id of the infected cell type according to a cc3d simulation
        """
        return getattr(self, self._infected_type_name.upper())

    def set_uninfected_type_name(self, _name: str):
        """
        Set the uninfected cell type name

        :param _name: Uninfected cell type name
        :return: None
        """
        self._uninfected_type_name = _name

    @property
    def uninfected_type_name(self) -> str:
        """
        Uninfected cell type name
        """
        return self._uninfected_type_name

    @uninfected_type_name.setter
    def uninfected_type_name(self, _name: str):
        self.set_uninfected_type_name(_name)

    def set_infected_type_name(self, _name: str):
        """
        Set the infected cell type name

        :param _name: Infected cell type name
        :return: None
        """
        self._infected_type_name = _name

    @property
    def infected_type_name(self) -> str:
        """
        Infected cell type name
        """
        return self._infected_type_name

    @infected_type_name.setter
    def infected_type_name(self, _name: str):
        self.set_infected_type_name(_name)

    def set_cell_diameter(self, _val: float):
        """
        Set the target cell diameter, in units of microns.
        Target volume for cells is interpreted as this value squared.

        :param _val: target cell diameter
        :return: None
        """
        if _val < 0:
            raise ValueError('Cell target diameter must be non-negative')
        self._cell_diameter = _val

    @property
    def cell_diameter(self) -> float:
        """
        Target cell diameter, in units of microns

        :return: None
        """
        return self._cell_diameter

    @cell_diameter.setter
    def cell_diameter(self, _val: float):
        self.set_cell_diameter(_val)

    @property
    def cell_volume(self):
        """
        Target volume for newly created cells; calculated from cell diameter
        """
        return self.cell_diameter * self.cell_diameter

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
    def volume_parameter(self, _val):
        self.set_volume_parameter(_val)


class VirusFieldInitializerSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Initializes virus field data and properties

    By default, requires CC3DML ids "virus_dc" for virus field diffusion coefficient and "virus_decay" for virus field
    decay. These can be set with set_field_data.
    """

    unique_key = virus_field_initializer_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self.virus_field_name = virus_field_name

        self._virus_diffusion_id = 'virus_dc'
        self._virus_decay_id = 'virus_decay'
        self._diffusion_coefficient = None
        self._decay_coefficient = None

    def start(self):
        """
        Called once to initialize simulation
        """
        diff_factor = self.step_period / (self.voxel_length * self.voxel_length)
        decay_factor = self.step_period
        if self._diffusion_coefficient is None:
            self.get_xml_element(self._virus_diffusion_id).cdata = ViralInfectionVTMModelInputs.virus_dc * diff_factor
        else:
            self.get_xml_element(self._virus_diffusion_id).cdata = self._diffusion_coefficient * diff_factor
        if self._decay_coefficient is None:
            self.get_xml_element(self._virus_decay_id).cdata = ViralInfectionVTMModelInputs.virus_decay * decay_factor
        else:
            self.get_xml_element(self._virus_decay_id).cdata = self._decay_coefficient * decay_factor

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
    def diffusion_coefficient(self) -> float:
        """
        Diffusion coefficient, in units microns^2/s
        """
        return self._diffusion_coefficient

    @diffusion_coefficient.setter
    def diffusion_coefficient(self, _val: float):
        if _val <= 0:
            raise ValueError("Diffusion coefficient must be positive")
        self._diffusion_coefficient = _val

    def set_decay_coefficient(self, _val: float):
        """
        Set decay coefficient, in units 1/s

        :param _val: Decay coefficient
        :return: None
        """
        if _val < 0.0:
            raise ValueError("Decay coefficient must be non-negative")
        self._decay_coefficient = _val

    @property
    def decay_coefficient(self) -> float:
        """
        Decay coefficient, in units 1/s
        """
        return self._decay_coefficient

    @decay_coefficient.setter
    def decay_coefficient(self, _val: float):
        self.set_diffusion_coefficient(_val)

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        """
        Set diffusion field data for virus field

        :param field_name: name of the virus field (optional)
        :param diffusion: cc3dml id of the diffusion field diffusion coefficient (optional)
        :param decay: cc3dml id of the diffusion field decay coefficient (optional)
        :return: None
        """
        if field_name is not None:
            self.virus_field_name = field_name
        if diffusion is not None:
            self._virus_diffusion_id = diffusion
        if decay is not None:
            self._virus_decay_id = decay

    @property
    def field_secretor(self):
        """
        Virus field secretor
        """
        return self.get_field_secretor(self.virus_field_name)

    @property
    def field_object(self):
        """
        Reference to virus field
        """
        return getattr(self.field, self.virus_field_name)


class ViralInternalizationSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral internalization

    The name of the infecting field can be set with the attribute target_field_name.

    Infection occurs for cells of the uninfected type, which converts them to the infected type.

    The names of the uninfected and infected types can be set with the attributes uninfected_type_name and
    infected_type_name, respectively.

    The rate of internalization can be set with the attribute internalization_rate
    """

    unique_key = viral_internalization_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self._target_field_name = ''
        self._uninfected_type_name = ''
        self._infected_type_name = ''
        self._internalization_rate = 0.0

        self.set_target_field_name(virus_field_name)
        self.set_uninfected_type_name(uninfected_type_name)
        self.set_infected_type_name(infected_type_name)
        self.set_internalization_rate(ViralInfectionVTMModelInputs.internalization_rate)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        secretor = self.get_field_secretor(field_name=self._target_field_name)
        for cell in self.cell_list_by_type(self.uninfected_type_id):
            seen_amount = secretor.amountSeenByCell(cell) / cell.volume
            rate = seen_amount * self._internalization_rate * self.step_period
            if random.random() <= nCoVUtils.ul_rate_to_prob(rate):
                self.set_cell_type(cell, self.infected_type_id)

    @property
    def uninfected_type_id(self) -> int:
        """
        Id of the uninfected cell type according to a cc3d simulation
        """
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def infected_type_id(self) -> int:
        """
        Id of the infected cell type according to a cc3d simulation
        """
        return getattr(self, self._infected_type_name.upper())

    def set_target_field_name(self, _name: str):
        """
        Set name of infecting field

        :param _name: Name of infecting field
        :return: None
        """
        self._target_field_name = _name

    @property
    def target_field_name(self):
        """
        Name of infecting field
        """
        return self._target_field_name

    @target_field_name.setter
    def target_field_name(self, _name: str):
        self.set_target_field_name(_name)

    def set_uninfected_type_name(self, _name: str):
        """
        Set the uninfected cell type name

        :param _name: Uninfected cell type name
        :return: None
        """
        self._uninfected_type_name = _name

    @property
    def uninfected_type_name(self):
        """
        Uninfected cell type name
        """
        return self._uninfected_type_name

    @uninfected_type_name.setter
    def uninfected_type_name(self, _name: str):
        self.set_uninfected_type_name(_name)

    def set_infected_type_name(self, _name: str):
        """
        Set the infected cell type name

        :param _name: Infected cell type name
        :return: None
        """
        self._infected_type_name = _name

    @property
    def infected_type_name(self):
        """
        Infected cell type name
        """
        return self._infected_type_name

    @infected_type_name.setter
    def infected_type_name(self, _name: str):
        self.set_infected_type_name(_name)

    def set_internalization_rate(self, _val: float):
        """
        Set internalization rate, in units virus/cell/s

        :param _val: Internalization rate
        :return: None
        """
        if _val < 0:
            raise ValueError("Internalization rate must be non-negative")
        self._internalization_rate = _val

    @property
    def internalization_rate(self):
        """
        Internalization rate, in units virus/cell/s
        """
        return self._internalization_rate

    @internalization_rate.setter
    def internalization_rate(self, _val):
        self.set_internalization_rate(_val)


class EclipsePhaseSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral eclipse phase

    After the eclipse phase, cells of the infected type are converted to the virus-releasing type.

    The names of the infected and virus-releasing types can be set with the attributes infected_type_name and
    virus_releasing_type_name, respectively.

    The period of the eclipse phase can be set with the attribute eclipse_phase.
    """

    unique_key = eclipse_phase_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self._infected_type_name = ''
        self._virus_releasing_type_name = ''
        self._eclipse_phase = 1.0

        self.set_infected_type_name(infected_type_name)
        self.set_virus_releasing_type_name(virus_releasing_type_name)
        self.set_eclipse_phase(ViralInfectionVTMModelInputs.eclipse_phase)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        pr = nCoVUtils.ul_rate_to_prob(self.step_period / self.eclipse_phase)
        for cell in self.cell_list_by_type(self.infected_type_id):
            if random.random() <= pr:
                self.set_cell_type(cell, self.virus_releasing_type_id)

    @property
    def infected_type_id(self) -> int:
        """
        Id of the infected cell type according to a cc3d simulation
        """
        return getattr(self, self._infected_type_name.upper())

    @property
    def virus_releasing_type_id(self) -> int:
        """
        Id of the virus-releasing cell type according to a cc3d simulation
        """
        return getattr(self, self._virus_releasing_type_name.upper())

    def set_infected_type_name(self, _name: str):
        """
        Set the infected cell type name

        :param _name: Infected cell type name
        :return: None
        """
        self._infected_type_name = _name

    @property
    def infected_type_name(self):
        """
        Infected cell type name
        """
        return self._infected_type_name

    @infected_type_name.setter
    def infected_type_name(self, _name: str):
        self.set_infected_type_name(_name)

    def set_virus_releasing_type_name(self, _name: str):
        """
        Set virus-releasing type name

        :param _name: Virus-releasing type name
        :return: None
        """
        self._virus_releasing_type_name = _name

    @property
    def virus_releasing_type_name(self):
        """
        Virus-releasing type name
        """
        return self._virus_releasing_type_name

    @virus_releasing_type_name.setter
    def virus_releasing_type_name(self, _name: str):
        self.set_virus_releasing_type_name(_name)

    def set_eclipse_phase(self, _val: float):
        """
        Set eclipse phase, in units of seconds

        :param _val: Eclipse phase
        :return: None
        """
        if _val <= 0:
            raise ValueError("Eclipse phase must be positive")
        self._eclipse_phase = _val

    @property
    def eclipse_phase(self):
        """
        Eclipse phase, in units of seconds
        """
        return self._eclipse_phase

    @eclipse_phase.setter
    def eclipse_phase(self, _val: float):
        self.set_eclipse_phase(_val)


class ViralDeathSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral death

    Cells of the virus-releasing type are converted to the dead type.

    The names of the virus-releasing and dead types can be set with the attributes virus_releasing_type_name and
    dead_type_name, respectively.

    The rate of viral death can be set with the attribute viral_death_rate.
    """

    unique_key = viral_death_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self._virus_releasing_type_name = ''
        self._dead_type_name = ''
        self._viral_death_rate = 0.0

        self.set_virus_releasing_type_name(virus_releasing_type_name)
        self.set_dead_type_name(dead_type_name)
        self.set_viral_death_rate(ViralInfectionVTMModelInputs.viral_death_rate)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        pr = nCoVUtils.ul_rate_to_prob(self.viral_death_rate * self.step_period)
        for cell in self.cell_list_by_type(self.virus_releasing_type_id):
            if random.random() <= pr:
                self.set_cell_type(cell, self.dead_type_id)

    @property
    def virus_releasing_type_id(self) -> int:
        """
        Id of the virus-releasing cell type according to a cc3d simulation
        """
        return getattr(self, self._virus_releasing_type_name.upper())

    @property
    def dead_type_id(self) -> int:
        """
        Id of the dead cell type according to a cc3d simulation
        """
        return getattr(self, self._dead_type_name.upper())

    def set_virus_releasing_type_name(self, _name: str):
        """
        Set the name of the virus-releasing cell type according to a cc3d simulation

        :param _name: Name of the virus-releasing cell type according to a cc3d simulation
        :return: None
        """
        self._virus_releasing_type_name = _name

    @property
    def virus_releasing_type_name(self):
        """
        Name of the virus-releasing cell type
        """
        return self._virus_releasing_type_name

    @virus_releasing_type_name.setter
    def virus_releasing_type_name(self, _name: str):
        self.set_virus_releasing_type_name(_name)

    def set_dead_type_name(self, _name: str):
        """
        Set the name of the dead cell type

        :param _name: Name of the dead cell type
        :return: None
        """
        self._dead_type_name = _name

    @property
    def dead_type_name(self):
        """
        Name of the dead cell type
        """
        return self._dead_type_name

    @dead_type_name.setter
    def dead_type_name(self, _name: str):
        self.set_dead_type_name(_name)

    def set_viral_death_rate(self, _val: float):
        """
        Set the rate of viral death, in units 1/s

        :param _val: Rate of viral death
        :return: None
        """
        if _val < 0:
            raise ValueError("Death rate must be non-negative")
        self._viral_death_rate = _val

    @property
    def viral_death_rate(self):
        """
        Rate of viral death, in units 1/s
        """
        return self._viral_death_rate

    @viral_death_rate.setter
    def viral_death_rate(self, _val: float):
        self.set_viral_death_rate(_val)


class ViralReleaseSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral release

    If the simulation domain is quasi-2D, then release will be multiplied by the height of the domain.

    The name of the virus-releasing cell type can be set with the attribute virus_releasing_type_name.

    The rate of virus release can be set with the attribute release_rate.

    The name of the field into which virus-releasing cells release virus can be set with the attribute
    target_field_name.
    """

    unique_key = viral_release_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self.runBeforeMCS = 1

        self._target_field_name = ''
        self._virus_releasing_type_name = ''
        self._release_rate = 0.0

        self.set_target_field_name(virus_field_name)
        self.set_virus_releasing_type_name(virus_releasing_type_name)
        self.set_release_rate(ViralInfectionVTMModelInputs.secretion_rate)

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        secretor = self.get_field_secretor(field_name=self._target_field_name)
        min_dim = min(self.dim.x, self.dim.y, self.dim.z)
        fact = 1.0
        if min_dim < 3:
            fact = float(min_dim)

        for cell in self.cell_list_by_type(self.virus_releasing_type_id):
            secretor.secreteInsideCell(cell, self.release_rate * fact / cell.volume * self.step_period)

    @property
    def virus_releasing_type_id(self) -> int:
        """
        Id of the virus-releasing cell type according to a cc3d simulation
        """
        return getattr(self, self._virus_releasing_type_name.upper())

    def set_target_field_name(self, _name: str):
        """
        Set the name of virus field

        :param _name: Name of virus field
        :return: None
        """
        self._target_field_name = _name

    @property
    def target_field_name(self):
        """
        Name of virus field
        """
        return self._target_field_name

    @target_field_name.setter
    def target_field_name(self, _name: str):
        self.set_target_field_name(_name)

    def set_virus_releasing_type_name(self, _name: str):
        """
        Set virus-releasing type name

        :param _name: Virus-releasing type name
        :return: None
        """
        self._virus_releasing_type_name = _name

    @property
    def virus_releasing_type_name(self):
        """
        Virus-releasing type name
        """
        return self._virus_releasing_type_name

    @virus_releasing_type_name.setter
    def virus_releasing_type_name(self, _name: str):
        self.set_virus_releasing_type_name(_name)

    def set_release_rate(self, _val: float):
        """
        Set the virus release rate, in units virus/cell/s

        :param _val: Virus release rate
        :return: None
        """
        if _val < 0:
            raise ValueError("Viral release must be non-negative")
        self._release_rate = _val

    @property
    def release_rate(self):
        """
        Virus release rate, in units virus/cell/s
        """
        return self._release_rate

    @release_rate.setter
    def release_rate(self, _val: float):
        self.set_release_rate(_val)


class SimDataSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Plots/writes simulation data of interest

    The frequency of plotting population data in Player can be set with the attribute plot_pop_data_freq.
    Plotting is disabled when plot_pop_data_freq is set to zero.

    The frequency of writing population data to file can be set with the attribute write_pop_data_freq.
    Data is written to the cc3d output directory with the name 'pop_data.dat'.
    Writing is disabled when write_pop_data_freq is set to zero.

    The frequency of plotting diffusion field data in Player can be set with the attribute plot_med_diff_data_freq.
    Plotting is disabled when plot_med_diff_data_freq is set to zero.

    The frequency of writing diffusion field data to file can be set with the attribute write_med_diff_data_freq.
    Data is written to the cc3d output directory with the name 'med_diff_data.dat'.
    Writing is disabled when write_med_diff_data_freq is set to zero.

    The name of each tracked cell type and the virus field can be set with the following attributes
        - uninfected: uninfected_type_name
        - infected: infected_type_name
        - virus releasing: virus_releasing_type_name
        - dead: dead_type_name
        - virus field: virus_field_name

    All data writing is performed every quarter-simulation period.
    """

    unique_key = sim_data_key

    def __init__(self, frequency=1):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self.pop_data_win = None
        self.pop_data_path = None
        self.pop_data = dict()

        self.med_diff_data_win = None
        self.med_diff_data_path = None
        self.med_diff_data = dict()

        self._plot_pop_data_freq = 0
        self._write_pop_data_freq = 0
        self._plot_med_diff_data_freq = 0
        self._write_med_diff_data_freq = 0

        self.plot_pop_data = False
        self.write_pop_data = False
        self.plot_med_diff_data = False
        self.write_med_diff_data = False
        self.med_diff_key = "MedDiff"

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

        self._uninfected_type_name = ''
        self._infected_type_name = ''
        self._virus_releasing_type_name = ''
        self._dead_type_name = ''
        self._virus_field_name = ''

        self.set_uninfected_type_name(uninfected_type_name)
        self.set_infected_type_name(infected_type_name)
        self.set_virus_releasing_type_name(virus_releasing_type_name)
        self.set_dead_type_name(dead_type_name)
        self.set_virus_field_name(virus_field_name)
        self.set_plot_pop_data_freq(ViralInfectionVTMModelInputs.plot_pop_data_freq)
        self.set_write_pop_data_freq(ViralInfectionVTMModelInputs.write_pop_data_freq)
        self.set_plot_med_diff_data_freq(ViralInfectionVTMModelInputs.plot_med_diff_data_freq)
        self.set_write_med_diff_data_freq(ViralInfectionVTMModelInputs.write_med_diff_data_freq)

    def start(self):
        """
        Called once to initialize simulation
        """
        self.plot_pop_data = self._plot_pop_data_freq > 0
        self.write_pop_data = self._write_pop_data_freq > 0

        self.plot_med_diff_data = self._plot_med_diff_data_freq > 0
        self.write_med_diff_data = self._write_med_diff_data_freq > 0

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

        if self.plot_med_diff_data:
            self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                              x_axis_title='MCS',
                                                              y_axis_title='Number of diffusive species per volume',
                                                              x_scale_type='linear',
                                                              y_scale_type='log',
                                                              grid=True,
                                                              config_options={'legend': True})

            self.med_diff_data_win.add_plot("MedViral", style='Dots', color='red', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath('pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

            if self.write_med_diff_data:
                self.med_diff_data_path = Path(self.output_dir).joinpath('med_diff_data.dat')
                with open(self.med_diff_data_path, 'w'):
                    pass

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """

        plot_pop_data = self.plot_pop_data and mcs % self._plot_pop_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % self._plot_med_diff_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % self._write_pop_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % self._write_med_diff_data_freq == 0
        else:
            write_pop_data = False
            write_med_diff_data = False

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.uninfected_type_id))
            num_cells_infected = len(self.cell_list_by_type(self.infected_type_id))
            num_cells_virusreleasing = len(self.cell_list_by_type(self.virus_releasing_type_id))
            num_cells_dying = len(self.cell_list_by_type(self.dead_type_id))

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

            # Write population data to file if requested
            if write_pop_data:
                self.pop_data[mcs] = [num_cells_uninfected,
                                      num_cells_infected,
                                      num_cells_virusreleasing,
                                      num_cells_dying]

        if plot_med_diff_data or write_med_diff_data:

            # Gather total diffusive amounts
            try:
                med_viral_total = self.get_field_secretor(self._virus_field_name).totalFieldIntegral()
            except AttributeError:  # Pre-v4.2.1 CC3D
                med_viral_total = 0.0
                field = getattr(self.field, self._virus_field_name)
                for x, y, z in self.every_pixel():
                    med_viral_total += field[x, y, z]

            # Plot total diffusive viral amount if requested
            if plot_med_diff_data:
                if med_viral_total > 0:
                    self.med_diff_data_win.add_data_point("MedViral", mcs, med_viral_total)

            # Write total diffusive viral amount if requested
            if write_med_diff_data:
                self.med_diff_data[mcs] = [med_viral_total]

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

    @property
    def uninfected_type_id(self) -> int:
        """
        Id of the uninfected cell type according to a cc3d simulation
        """
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def infected_type_id(self) -> int:
        """
        Id of the infected cell type according to a cc3d simulation
        """
        return getattr(self, self._infected_type_name.upper())

    @property
    def virus_releasing_type_id(self) -> int:
        """
        Id of the virus-releasing cell type according to a cc3d simulation
        """
        return getattr(self, self._virus_releasing_type_name.upper())

    @property
    def dead_type_id(self) -> int:
        """
        Id of the dead cell type according to a cc3d simulation
        """
        return getattr(self, self._dead_type_name.upper())

    def set_uninfected_type_name(self, _name: str):
        """
        Set the uninfected cell type name

        :param _name: Uninfected cell type name
        :return: None
        """
        self._uninfected_type_name = _name

    @property
    def uninfected_type_name(self):
        """
        Uninfected cell type name
        """
        return self._uninfected_type_name

    @uninfected_type_name.setter
    def uninfected_type_name(self, _name: str):
        self.set_uninfected_type_name(_name)

    def set_infected_type_name(self, _name: str):
        """
        Set the infected cell type name

        :param _name: Infected cell type name
        :return: None
        """
        self._infected_type_name = _name

    @property
    def infected_type_name(self):
        """
        Infected cell type name
        """
        return self._infected_type_name

    @infected_type_name.setter
    def infected_type_name(self, _name: str):
        self.set_infected_type_name(_name)

    def set_virus_releasing_type_name(self, _name: str):
        """
        Set virus-releasing type name

        :param _name: Virus-releasing type name
        :return: None
        """
        self._virus_releasing_type_name = _name

    @property
    def virus_releasing_type_name(self):
        """
        Virus-releasing type name
        """
        return self._virus_releasing_type_name

    @virus_releasing_type_name.setter
    def virus_releasing_type_name(self, _name: str):
        self.set_virus_releasing_type_name(_name)

    def set_dead_type_name(self, _name: str):
        """
        Set the name of the dead cell type

        :param _name: Name of the dead cell type
        :return: None
        """
        self._dead_type_name = _name

    @property
    def dead_type_name(self):
        """
        Name of the dead cell type
        """
        return self._dead_type_name

    @dead_type_name.setter
    def dead_type_name(self, _name: str):
        self.set_dead_type_name(_name)

    def set_virus_field_name(self, _name: str):
        """
        Set the name of the virus field

        :param _name: Name of the virus field
        :return: None
        """
        self._virus_field_name = _name

    @property
    def virus_field_name(self):
        """
        Name of the virus field
        """
        return self._virus_field_name

    @virus_field_name.setter
    def virus_field_name(self, _name: str):
        self.set_virus_field_name(_name)

    def set_plot_pop_data_freq(self, _val: int):
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

    def set_write_med_diff_data_freq(self, _val: int):
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
                       (self.write_med_diff_data, self.med_diff_data_path, self.med_diff_data)]
        for write_data, data_path, data in output_info:
            if write_data:
                with open(data_path, 'a') as fout:
                    fout.write(self.data_output_string(data))
                    data.clear()
