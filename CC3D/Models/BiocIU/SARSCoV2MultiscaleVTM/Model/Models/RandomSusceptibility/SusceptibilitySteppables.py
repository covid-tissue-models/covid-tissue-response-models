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

from random import random

from ViralInfectionVTMSteppables import uninfected_type_name
from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
from Models.SegoAponte2020.ViralInfectionVTMLib import unbound_receptors_cellg_key

from .SusceptibilityModelInputs import *

vs_state_key = 'vs_susc'
random_susceptibility_steppable_key = 'random_susceptibility_steppable'  # RandomSusceptibilitySteppable


class RandomSusceptibilitySteppable(nCoVSteppableBase):
    """
    Implements randomly varying viral susceptibility in space
    Implementation of non-susceptibility can be overwritten in subclasses (default sets surface receptors to zero)
    Random susceptibility is implemented by assigning a probability of a cell being initialized as unsusceptible
    according to a fractional value of unsusceptible cells.
    """

    unique_key = random_susceptibility_steppable_key

    def __init__(self, frequency=1):
        nCoVSteppableBase.__init__(self, frequency)

        if track_susc:
            self.track_cell_level_scalar_attribute(field_name='Susceptible', attribute_name=vs_state_key)

        self._epithelial_types = {}

        # Initialize default data
        self.append_epithelial_type(uninfected_type_name, frac_not_susc)

    def make_unsusceptible(self, _cell) -> bool:
        """
        Try to make cell not susceptible to infection; can be overwritten in subclasses
        :param _cell: cell to try to make not susceptible
        :return: True if cell is made not susceptible; False otherwise
        """
        assert unbound_receptors_cellg_key in _cell.dict.keys(), \
            'RandomSusceptibilitySteppable requires surface receptors'
        if _cell.dict[unbound_receptors_cellg_key] > 0:
            _cell.dict[unbound_receptors_cellg_key] = 0
            if track_susc:
                _cell.dict[vs_state_key] = False
            return True
        return False

    def append_epithelial_type(self, _name: str, _val: float = frac_not_susc):
        if not 0 <= _val <= 1:
            raise ValueError('Fractional values must be in [0, 1]')
        self._epithelial_types[_name] = _val

    def remove_epithelial_type(self, _name: str):
        self._epithelial_types.pop(_name)

    @property
    def epithelial_type_ids(self):
        return [getattr(self, x.upper()) for x in self._epithelial_types.keys()]

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback. If a cell's type is registered, then it becomes unsusceptible with a probability
        according to its registered fraction
        """
        if _new_cell.type in self.epithelial_type_ids:
            if track_susc:
                _new_cell.dict[vs_state_key] = True
            if random() <= self._epithelial_types[self.get_type_name_by_cell(_new_cell)]:
                self.make_unsusceptible(_new_cell)
