"""
Defines module steppables

Steppables
==========

SimpleRecoverySteppable
-----------------------
Description: implements recovery model

Usage: In ViralInfectionVTM.py, add the following

from Models.RecoverySimple.RecoverySteppables import SimpleRecoverySteppable

CompuCellSetup.register_steppable(steppable=SimpleRecoverySteppable(frequency=1))

SimpleRecoveryDataSteppable
---------------------------

Description: performs data tracking

Usage: In ViralInfectionVTM.py, add the following

from Models.RecoverySimple.RecoverySteppables import SimpleRecoveryDataSteppable

CompuCellSetup.register_steppable(steppable=SimpleRecoveryDataSteppable(frequency=1))
"""

from collections import namedtuple
import random
import os
import sys
from typing import *

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
from ViralInfectionVTMSteppables import SimDataSteppable

from . import RecoveryInputs
rec_steppable_key = "sprec_steppable"
rec_data_steppable_key = "sprec_data_steppable"


RecoveryData = namedtuple(typename='RecoveryData', field_names=['recovered_type_name', 'recovery_rate'])


class SimpleRecoverySteppable(nCoVSteppableBase):
    """
    Implements simple recovery
    """

    unique_key = rec_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.num_recovered = 0
        self.num_recovered_by_type: Dict[str, int] = {}

        self._registered_recovery: Dict[str, RecoveryData] = {}

        # Initialize default data
        from ViralInfectionVTMSteppables import uninfected_type_name, dead_type_name
        self.register_recovery(dead_type_name=dead_type_name,
                               recovered_type_name=uninfected_type_name,
                               recovery_rate=RecoveryInputs.recovery_rate)

    def step(self, mcs):
        [self.recover_cell(cell) for cell in self.cell_list_by_type(*self.dead_type_ids) if self.cell_recovers(cell)]

    def cell_recovers(self, _cell) -> bool:
        """
        Recovery criterion

        :param _cell: dead cell to test for recovery
        :return: True if cell recovers
        """
        return random.random() < self.recovery_prob(_cell)

    def recover_cell(self, _cell):
        """
        Implement recovery

        :param _cell: dead cell to recover
        :return: None
        """
        self.num_recovered += 1
        self.num_recovered_by_type[self.get_type_name_by_cell(_cell)] += 1
        self.set_cell_type(_cell, self.recovery_type(_cell))

    @property
    def dead_type_ids(self) -> List[int]:
        """
        Ids of dead cell types according to a cc3d simulation
        """
        return [getattr(self, x.upper()) for x in self._registered_recovery.keys()]

    def recovery_prob(self, _cell) -> float:
        """
        Returns probability of recovery for a cell

        :param _cell: a cell
        :return: Probability of recovery for a cell
        """
        return self._registered_recovery[self.get_type_name_by_cell(_cell)].recovery_rate * self.step_period

    def recovery_type(self, _cell) -> int:
        """
        Id of recovered type that corresponds to a cell according to a cc3d simulation

        :param _cell: a cell
        :return: Id of recovered type
        """
        return getattr(self, self._registered_recovery[self.get_type_name_by_cell(_cell)].recovered_type_name.upper())

    def register_recovery(self, dead_type_name: str, recovered_type_name: str, recovery_rate: float):
        """
        Register a recovery data entry

        :param dead_type_name: name of dead cell type
        :param recovered_type_name: name of recovered cell type
        :param recovery_rate: rate of recovery (1/s)
        :return: None
        """
        self._registered_recovery[dead_type_name] = RecoveryData(recovered_type_name=recovered_type_name,
                                                                 recovery_rate=recovery_rate)
        self.num_recovered_by_type[dead_type_name] = 0

    def unregister_recovery(self, _name: str):
        """
        Unregister a recovery data entry

        :param _name: name of dead cell type
        :return: None
        """
        self._registered_recovery.pop(_name)
        self.num_recovered_by_type.pop(_name)


class SimpleRecoveryDataSteppable(nCoVSteppableBase):
    """
    Implements simple recovery data tracking; like SimDataSteppable in main framework
    """

    unique_key = rec_data_steppable_key

    def __init__(self, frequency=1, plot_freq=None, write_freq=None):
        super().__init__(frequency)

        self.rec_steppable = None

        self.rec_data_win = None
        self.rec_data_path = None
        self.rec_data = dict()

        self._flush_counter = 1

        self.plot_rec_data_freq = RecoveryInputs.plot_rec_data_freq
        self.write_rec_data_freq = RecoveryInputs.write_rec_data_freq

        # Subclass support
        if plot_freq is not None:
            self.plot_rec_data_freq = plot_freq
        if write_freq is not None:
            self.write_rec_data_freq = write_freq

        # Particularization; subclasses should set these to something unique
        self.rec_steppable_key = rec_steppable_key
        self.data_path_rel = 'sprec_data.dat'
        self.plot_title = 'Simple Recovery'

    def start(self):
        """
        Called once to initialize simulation
        """
        # Initialze plot window if requested
        if self.plot_rec_data_freq > 0:
            self.rec_data_win = self.add_new_plot_window(title=self.plot_title,
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Recovered cells', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=False,
                                                         config_options={'legend': True})

            self.rec_data_win.add_plot("Recovered", style='Dots', color='blue', size=5)

        # Initialize data output if requested
        from pathlib import Path
        if self.write_rec_data_freq > 0:
            self.rec_data_path = Path(self.output_dir).joinpath(self.data_path_rel)
            with open(self.rec_data_path, 'w'):
                pass

    def step(self, mcs):
        """
        Called every simulation step

        :param mcs: current simulation step
        :return: None
        """
        if self.rec_steppable is None:
            self.rec_steppable = self.shared_steppable_vars[self.rec_steppable_key]

        # Plot population data plot if requested
        if self.plot_rec_data_freq > 0 and mcs % self.plot_rec_data_freq == 0:
            self.rec_data_win.add_data_point("Recovered", mcs, self.rec_steppable.num_recovered)

        # Write population data to file if requested
        if self.write_rec_data_freq > 0 and mcs % self.write_rec_data_freq == 0:
            self.rec_data[mcs] = [self.rec_steppable.num_recovered]

        # Flush outputs at quarter simulation lengths
        if mcs >= int(self.simulator.getNumSteps() / 4 * self._flush_counter):
            self.flush_stored_outputs()
            self._flush_counter += 1

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

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage

        :return: None
        """
        if self.write_rec_data_freq > 0:
            with open(self.rec_data_path, 'a') as fout:
                fout.write(SimDataSteppable.data_output_string(self, self.rec_data))
                self.rec_data.clear()
