# Model of neighbor-dependent recovery
# Written by T.J. Sego, Ph.D., and presented in the Interactive Two-Part Virtual Mini-workshop on
# Open-Source CompuCell3D Multiscale, Virtual-Tissue Spatiotemporal Modeling and Simulations of COVID-19 Infection,
# Viral Spread and Immune Response and Treatment Regimes, June 11-12 & Jun 18-19, 2020.
# Dead cells "resurrect" and become uninfected according to the number of uninfected cells in its neighborhood and a
# rate "recovery_rate" defined in RecoveryInputs.py.
#
# Built from RecoverySimple module written by T.J. Sego, Ph.D.
#
# NeighborRecoverySteppable
#   Description: implements recovery model
#   Usage:
#       In ViralInfectionVTM.py, add the following
#           from Models.RecoveryNeighbor.RecoverySteppables import NeighborRecoverySteppable
#           CompuCellSetup.register_steppable(steppable=NeighborRecoverySteppable(frequency=1))
# NeighborRecoveryDataSteppable
#   Description: performs data tracking
#   Usage:
#       In ViralInfectionVTM.py, add the following
#           from Models.RecoveryNeighbor.RecoverySteppables import NeighborRecoveryDataSteppable
#           CompuCellSetup.register_steppable(steppable=NeighborRecoveryDataSteppable(frequency=1))

import random
import sys
import os
from cc3d.core.PySteppables import *

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
from ViralInfectionVTMModelInputs import s_to_mcs

from .RecoveryInputs import *
rec_steppable_key = "nbrec_steppable"

# Inherit from Simple Recovery model
sys.path.append(os.environ["ViralInfectionVTM"])
from Models.RecoverySimple.RecoverySteppables import SimpleRecoverySteppable, SimpleRecoveryDataSteppable


class NeighborRecoverySteppable(SimpleRecoverySteppable):
    """
    Implements neighbor-dependent recovery
    """
    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Particularize SimpleRecoverySteppable
        self.rec_steppable_key = rec_steppable_key

    def cell_recovers(self, _cell) -> bool:
        """
        Recovery criterion
        :param _cell: dead cell to test for recovery
        :return: True if cell recovers
        """
        ca = sum([a for n, a in self.get_cell_neighbor_data_list(_cell) if n is not None and n.type == self.UNINFECTED])
        return random.random() < ca * recovery_rate * s_to_mcs


class NeighborRecoveryDataSteppable(SimpleRecoveryDataSteppable):
    """
    Implements neighbor-dependent recovery data tracking; like SimDataSteppable in main framework
    """
    def __init__(self, frequency=1):
        super().__init__(frequency, plot_freq=plot_rec_data_freq, write_freq=write_rec_data_freq)

        # Particularize SimpleRecoveryDataSteppable
        self.rec_steppable_key = rec_steppable_key
        self.data_path_rel = 'nbrec_data.dat'
        self.plot_title = 'Neighbor Recovery'
