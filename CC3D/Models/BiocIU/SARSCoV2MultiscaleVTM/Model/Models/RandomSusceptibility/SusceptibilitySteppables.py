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

from random import randint

from ViralInfectionVTMSteppables import (
    uninfected_type_name, infected_type_name, virus_releasing_type_name, dead_type_name)
from Models.SegoAponte2020 import ViralInfectionVTMBasePy
from Models.SegoAponte2020.ViralInfectionVTMLib import unbound_receptors_cellg_key

from .SusceptibilityModelInputs import *

ViralInfectionVTMSteppableBasePy = ViralInfectionVTMBasePy.ViralInfectionVTMSteppableBasePy

vs_state_key = 'vs_susc'
random_susceptibility_steppable_key = 'random_susceptibility_steppable'  # RandomSusceptibilitySteppable


class RandomSusceptibilitySteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements randomly varying viral susceptibility in space
    Implementation of non-susceptibility can be overwritten in subclasses (default sets surface receptors to zero)
    """

    unique_key = random_susceptibility_steppable_key

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        if track_susc:
            self.track_cell_level_scalar_attribute(field_name='Susceptible', attribute_name=vs_state_key)

        self._epithelial_types = []
        self._target_type = ''
        self._frac_not_susc = 0.0

        # Initialize default data
        self.append_epithelial_type(uninfected_type_name)
        self.append_epithelial_type(infected_type_name)
        self.append_epithelial_type(virus_releasing_type_name)
        self.append_epithelial_type(dead_type_name)
        self.set_target_type(uninfected_type_name)
        self.set_frac_not_susc(frac_not_susc)

    def start(self):
        ec_list = self.cell_list_by_type(*[getattr(self, x.upper()) for x in self._epithelial_types])
        if track_susc:
            for cell in ec_list:
                cell.dict[vs_state_key] = True
        num_changed = 0
        while num_changed < int(len(ec_list) * self._frac_not_susc):
            if self.make_unsusceptible(self.cell_field[randint(0, self.dim.x - 1), randint(0, self.dim.y - 1), 0]):
                num_changed += 1

    def make_unsusceptible(self, _cell) -> bool:
        """
        Try to make cell not susceptible to infection; can be overwritten in subclasses
        :param _cell: cell to try to make not susceptible
        :return: True if cell is made not susceptible; False otherwise
        """
        assert unbound_receptors_cellg_key in _cell.dict.keys(), \
            'RandomSusceptibilitySteppable requires surface receptors'
        if _cell.dict[unbound_receptors_cellg_key] > 0 and _cell.type == self.target_type_id:
            _cell.dict[unbound_receptors_cellg_key] = 0
            if track_susc:
                _cell.dict[vs_state_key] = False
            return True
        return False

    def append_epithelial_type(self, _name: str):
        self._epithelial_types.append(_name)

    def remove_epithelial_type(self, _name: str):
        self._epithelial_types.remove(_name)

    def set_target_type(self, _name: str):
        self._target_type = _name

    @property
    def target_type_id(self) -> int:
        return getattr(self, self._target_type.upper())

    def set_frac_not_susc(self, _frac: float):
        assert 0 <= _frac < 1, 'Fraction of not susceptible cells (num_not_susc) must be in [0, 1)'
        self._frac_not_susc = _frac
