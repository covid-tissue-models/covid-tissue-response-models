###############################################################################################################
# To cite this model please use the following:
#
# T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
# Jessica Crawshaw, James M. Osborne, Ellen M. Quardokus, Richard K. Plemper, James A. Glazier,
# "A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and
# immune response in epithelial tissues and its application to drug therapy timing and effectiveness",
# PLoS Comput Biol 16(12): e1008451. https://doi.org/10.1371/journal.pcbi.1008451
###############################################################################################################

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
    Initializes an epithelial sheet
    """

    unique_key = cell_initializer_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self._initialization_mode = None
        self._infection_mode = None
        self._uninfected_type_name = ''
        self._infected_type_name = ''

        self._aux_infect_vars = None

        self.set_uninfected_type_name(uninfected_type_name)
        self.set_infected_type_name(infected_type_name)

        self.initialize_sheet()
        self.random_infected_fraction(0.01)

    def start(self):
        # Do cell initialization
        self.initialize_cells()

        # Do initial infection
        self.initialize_infection()

    def initialize_cells(self):
        if self._initialization_mode is not None:
            self._initialization_mode()

    def initialize_infection(self):
        if self._infection_mode is not None:
            self._infection_mode()

    def _initialize_sheet(self):
        # Enforce compatible lattice dimensions with epithelial cell size
        cdiam = ViralInfectionVTMModelInputs.cell_diameter
        assert self.dim.x % cdiam == 0 and self.dim.y % cdiam == 0, \
            f'Lattice dimensions must be multiples of the unitless cell diameter (currently cell_diameter = {cdiam})'

        for x in range(0, self.dim.x, int(cdiam)):
            for y in range(0, self.dim.y, int(cdiam)):
                cell = self.new_cell(self.uninfected_type_id)
                self.cellField[x:x + int(cdiam), y:y + int(cdiam), 0] = cell

                cell.targetVolume = ViralInfectionVTMModelInputs.cell_volume
                cell.lambdaVolume = ViralInfectionVTMModelInputs.volume_lm

    def initialize_sheet(self):
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
        self._infection_mode = self._no_initial_infection
        self._aux_infect_vars = None

    def single_infected_cell(self):
        self._infection_mode = self._single_infected_cell
        self._aux_infect_vars = None

    def random_infected_fraction(self, _frac: float):
        if not (0 <= _frac <= 1):
            raise ValueError("Invalid initial infection fraction: must be in [0, 1]")
        self._infection_mode = self._random_infected_fraction
        self._aux_infect_vars = _frac

    def random_infected_probability(self, _prob: float):
        if not (0 <= _prob <= 1):
            raise ValueError("Invalid initial infection probability: must be in [0, 1]")
        self._infection_mode = self._random_infected_probability
        self._aux_infect_vars = _prob

    @property
    def uninfected_type_id(self) -> int:
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def infected_type_id(self) -> int:
        return getattr(self, self._infected_type_name.upper())

    def set_uninfected_type_name(self, _name: str):
        self._uninfected_type_name = _name

    def set_infected_type_name(self, _name: str):
        self._infected_type_name = _name


class VirusFieldInitializerSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Initializes field data and properties

    By default, requires CC3DML ids "virus_dc" for virus field diffusion coefficient and "virus_decay" for virus field
    decay
    """

    unique_key = virus_field_initializer_key

    def __init__(self, frequency):
        nCoVSteppableBase.nCoVSteppableBase.__init__(self, frequency)

        self.virus_field_name = virus_field_name

        self._virus_diffusion_id = 'virus_dc'
        self._virus_decay_id = 'virus_decay'

    def start(self):
        self.diffusion_coefficient = ViralInfectionVTMModelInputs.virus_dc
        self.decay_coefficient = ViralInfectionVTMModelInputs.virus_decay

    @property
    def diffusion_coefficient(self) -> float:
        return self.get_xml_element(self._virus_diffusion_id).cdata

    @diffusion_coefficient.setter
    def diffusion_coefficient(self, _val: float):
        if _val <= 0:
            raise ValueError("Diffusion coefficient must be positive")
        self.get_xml_element(self._virus_diffusion_id).cdata = _val

    @property
    def decay_coefficient(self) -> float:
        return self.get_xml_element(self._virus_decay_id).cdata

    @decay_coefficient.setter
    def decay_coefficient(self, _val: float):
        if _val < 0:
            raise ValueError("Decay coefficient must be non-negative")
        self.get_xml_element(self._virus_decay_id).cdata = _val

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        if field_name is not None:
            self.virus_field_name = field_name
        if diffusion is not None:
            self._virus_diffusion_id = diffusion
        if decay is not None:
            self._virus_decay_id = decay

    @property
    def field_secretor(self):
        return self.get_field_secretor(self.virus_field_name)

    @property
    def field_object(self):
        return getattr(self.field, self.virus_field_name)


class ViralInternalizationSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral internalization
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
        secretor = self.get_field_secretor(field_name=self._target_field_name)
        for cell in self.cell_list_by_type(self.uninfected_type_id):
            seen_amount = secretor.amountSeenByCell(cell) / cell.volume
            rate = seen_amount * self._internalization_rate
            if random.random() <= nCoVUtils.ul_rate_to_prob(rate):
                self.set_cell_type(cell, self.infected_type_id)

    @property
    def uninfected_type_id(self) -> int:
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def infected_type_id(self) -> int:
        return getattr(self, self._infected_type_name.upper())

    def set_target_field_name(self, _name: str):
        self._target_field_name = _name

    def set_uninfected_type_name(self, _name: str):
        self._uninfected_type_name = _name

    def set_infected_type_name(self, _name: str):
        self._infected_type_name = _name

    def set_internalization_rate(self, _val: float):
        if _val < 0:
            raise ValueError("Internalization rate must be non-negative")
        self._internalization_rate = _val


class EclipsePhaseSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral eclipse phase
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
        pr = nCoVUtils.ul_rate_to_prob(1.0 / self._eclipse_phase)
        for cell in self.cell_list_by_type(self.infected_type_id):
            if random.random() <= pr:
                self.set_cell_type(cell, self.virus_releasing_type_id)

    @property
    def infected_type_id(self) -> int:
        return getattr(self, self._infected_type_name.upper())

    @property
    def virus_releasing_type_id(self) -> int:
        return getattr(self, self._virus_releasing_type_name.upper())

    def set_infected_type_name(self, _name: str):
        self._infected_type_name = _name

    def set_virus_releasing_type_name(self, _name: str):
        self._virus_releasing_type_name = _name

    def set_eclipse_phase(self, _val: float):
        if _val <= 0:
            raise ValueError("Eclipse phase must be positive")
        self._eclipse_phase = _val


class ViralDeathSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral death
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
        pr = nCoVUtils.ul_rate_to_prob(self._viral_death_rate)
        for cell in self.cell_list_by_type(self.virus_releasing_type_id):
            if random.random() <= pr:
                self.set_cell_type(cell, self.dead_type_id)

    @property
    def virus_releasing_type_id(self) -> int:
        return getattr(self, self._virus_releasing_type_name.upper())

    @property
    def dead_type_id(self) -> int:
        return getattr(self, self._dead_type_name.upper())

    def set_virus_releasing_type_name(self, _name: str):
        self._virus_releasing_type_name = _name

    def set_dead_type_name(self, _name: str):
        self._dead_type_name = _name

    def set_viral_death_rate(self, _val: float):
        if _val < 0:
            raise ValueError("Death rate must be non-negative")
        self._viral_death_rate = _val


class ViralReleaseSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Performs viral release

    If the simulation domain is quasi-2D, then release will be multiplied by the height of the domain
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
        secretor = self.get_field_secretor(field_name=self._target_field_name)
        min_dim = min(self.dim.x, self.dim.y, self.dim.z)
        fact = 1.0
        if min_dim < 3:
            fact = float(min_dim)

        for cell in self.cell_list_by_type(self.virus_releasing_type_id):
            secretor.secreteInsideCell(cell, self._release_rate * fact / cell.volume)

    @property
    def virus_releasing_type_id(self) -> int:
        return getattr(self, self._virus_releasing_type_name.upper())

    def set_target_field_name(self, _name: str):
        self._target_field_name = _name

    def set_virus_releasing_type_name(self, _name: str):
        self._virus_releasing_type_name = _name

    def set_release_rate(self, _val: float):
        if _val < 0:
            raise ValueError("Viral release must be non-negative")
        self._release_rate = _val


class SimDataSteppable(nCoVSteppableBase.nCoVSteppableBase):
    """
    Plots/writes simulation data of interest
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

        self.plot_pop_data = ViralInfectionVTMModelInputs.plot_pop_data_freq > 0
        self.write_pop_data = ViralInfectionVTMModelInputs.write_pop_data_freq > 0

        self.plot_med_diff_data = ViralInfectionVTMModelInputs.plot_med_diff_data_freq > 0
        self.write_med_diff_data = ViralInfectionVTMModelInputs.write_med_diff_data_freq > 0
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

    def start(self):
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

        plot_pop_data = self.plot_pop_data and mcs % ViralInfectionVTMModelInputs.plot_pop_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % ViralInfectionVTMModelInputs.plot_med_diff_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % ViralInfectionVTMModelInputs.write_pop_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % ViralInfectionVTMModelInputs.write_med_diff_data_freq == 0
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
        self.finish()

    def finish(self):
        self.flush_stored_outputs()

    @property
    def uninfected_type_id(self) -> int:
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def infected_type_id(self) -> int:
        return getattr(self, self._infected_type_name.upper())

    @property
    def virus_releasing_type_id(self) -> int:
        return getattr(self, self._virus_releasing_type_name.upper())

    @property
    def dead_type_id(self) -> int:
        return getattr(self, self._dead_type_name.upper())

    def set_uninfected_type_name(self, _name: str):
        self._uninfected_type_name = _name

    def set_infected_type_name(self, _name: str):
        self._infected_type_name = _name

    def set_virus_releasing_type_name(self, _name: str):
        self._virus_releasing_type_name = _name

    def set_dead_type_name(self, _name: str):
        self._dead_type_name = _name

    def set_virus_field_name(self, _name: str):
        self._virus_field_name = _name

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
