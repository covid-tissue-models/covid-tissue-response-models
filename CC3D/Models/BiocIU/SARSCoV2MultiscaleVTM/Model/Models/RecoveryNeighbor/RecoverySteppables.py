"""
Defines module steppables

Steppables
==========

NeighborRecoverySteppable
-------------------------
Description: implements recovery model

Usage: In ViralInfectionVTM.py, add the following

from Models.RecoveryNeighbor.RecoverySteppables import NeighborRecoverySteppable

CompuCellSetup.register_steppable(steppable=NeighborRecoverySteppable(frequency=1))

NeighborRecoveryDataSteppable
-----------------------------
Description: performs data tracking

Usage: In ViralInfectionVTM.py, add the following

from Models.RecoveryNeighbor.RecoverySteppables import NeighborRecoveryDataSteppable

CompuCellSetup.register_steppable(steppable=NeighborRecoveryDataSteppable(frequency=1))
"""

import random

from . import RecoveryInputs
rec_steppable_key = "nbrec_steppable"
rec_data_steppable_key = "nbrec_data_steppable"

# Inherit from Simple Recovery model
from Models.RecoverySimple.RecoverySteppables import SimpleRecoverySteppable, SimpleRecoveryDataSteppable


class NeighborRecoverySteppable(SimpleRecoverySteppable):
    """
    Implements neighbor-dependent recovery
    """

    unique_key = rec_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self._voxel_area = self.voxel_length * self.voxel_length

        # Initialize default data
        from ViralInfectionVTMSteppables import uninfected_type_name, dead_type_name
        self.register_recovery(dead_type_name=dead_type_name,
                               recovered_type_name=uninfected_type_name,
                               recovery_rate=RecoveryInputs.recovery_rate)

    def start(self):
        super().start()

        # Reducing a little computational cost here
        self._voxel_area = self.voxel_length * self.voxel_length

    def cell_recovers(self, _cell) -> bool:
        """
        Recovery criterion

        :param _cell: dead cell to test for recovery
        :return: True if cell recovers
        """
        recovery_type = self.recovery_type(_cell)
        ca = sum([a for n, a in self.get_cell_neighbor_data_list(_cell) if n is not None and n.type == recovery_type])
        if ca == 0:
            return False
        return random.random() < ca * self._voxel_area * self.recovery_prob(_cell)


class NeighborRecoveryDataSteppable(SimpleRecoveryDataSteppable):
    """
    Implements neighbor-dependent recovery data tracking; like SimDataSteppable in main framework
    """

    unique_keys = rec_data_steppable_key

    def __init__(self, frequency=1, plot_rec_data_freq=None, write_rec_data_freq=None):
        if plot_rec_data_freq is None:
            plot_rec_data_freq = RecoveryInputs.plot_rec_data_freq
        if write_rec_data_freq is None:
            write_rec_data_freq = RecoveryInputs.write_rec_data_freq
        super().__init__(frequency, plot_freq=plot_rec_data_freq, write_freq=write_rec_data_freq)

        # Particularize SimpleRecoveryDataSteppable
        self.rec_steppable_key = rec_steppable_key
        self.data_path_rel = 'nbrec_data.dat'
        self.plot_title = 'Neighbor Recovery'
