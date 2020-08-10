# Randomly distributed virus susceptibility
# Written by T.J. Sego, Ph.D.
# A fraction of uninfected cells are randomly selected and made unsusceptible to viral internalization
# Model parameters are specified in SusceptibilityModelInputs.py
#
# RandomSusceptibilitySteppable
#   Description: implements randomly distributed susceptibility
#   Usage:
#       In ViralInfectionVTM.py, add the following
#           from Models.RandomSusceptibility.SusceptibilitySteppables import RandomSusceptibilitySteppable
#           CompuCellSetup.register_steppable(steppable=RandomSusceptibilitySteppable(frequency=1))

import os
from random import randint
import sys
sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))

from ViralInfectionVTMSteppableBasePy import ViralInfectionVTMSteppableBasePy

from .SusceptibilityModelInputs import *

vs_state_key = 'vs_susc'

class RandomSusceptibilitySteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements randomly varying viral susceptibility in space
    Implementation of non-susceptibility can be overwritten in subclasses (default sets surface receptors to zero)
    """
    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        assert 0 <= frac_not_susc < 1, 'Fraction of not susceptible cells (num_not_susc) must be in [0, 1)'
        if track_susc:
            self.track_cell_level_scalar_attribute(field_name='Susceptible', attribute_name=vs_state_key)

    def start(self):
        ec_list = self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.VIRUSRELEASING, self.DYING)
        if track_susc:
            for cell in ec_list:
                cell.dict[vs_state_key] = True
        num_changed = 0
        while num_changed < int(len(ec_list) * frac_not_susc):
            if self.make_unsusceptible(self.cell_field[randint(0, self.dim.x - 1), randint(0, self.dim.y - 1), 0]):
                num_changed += 1

    def make_unsusceptible(self, _cell) -> bool:
        """
        Try to make cell not susceptible to infection; can be overwritten in subclasses
        :param _cell: cell to try to make not susceptible
        :return: True if cell is made not susceptible; False otherwise
        """
        assert 'Receptors' in _cell.dict.keys(), 'RandomSusceptibilitySteppable requires surface receptors'
        if _cell.dict['Receptors'] > 0 and _cell.type == self.UNINFECTED:
            _cell.dict['Receptors'] = 0
            if track_susc:
                _cell.dict[vs_state_key] = False
            return True
        return False
