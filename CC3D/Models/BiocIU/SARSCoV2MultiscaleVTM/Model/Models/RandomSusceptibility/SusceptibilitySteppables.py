"""
Defines module steppables

Steppables
==========

RandomSusceptibilitySteppable
-----------------------------
Description: implements randomly distributed susceptibility

Usage: In ViralInfectionVTM.py, add the following

from Models.RandomSusceptibility.SusceptibilitySteppables import RandomSusceptibilitySteppable

steppable = RandomSusceptibilitySteppable(frequency=1)

# Add a cell type called 'MyType' and assign a fraction of unsusceptible cells of 1/10.

steppable.append_epithelial_type('MyType', 0.01)

CompuCellSetup.register_steppable(steppable=steppable)
"""

from random import random

from ViralInfectionVTMSteppables import uninfected_type_name
from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
from Models.SegoAponte2020.ViralInfectionVTMLib import unbound_receptors_cellg_key

from . import SusceptibilityModelInputs

vs_state_key = 'vs_susc'
random_susceptibility_steppable_key = 'random_susceptibility_steppable'  # RandomSusceptibilitySteppable


class RandomSusceptibilitySteppable(nCoVSteppableBase):
    """
    Implements randomly varying viral susceptibility in space

    Implementation of non-susceptibility can be overwritten in subclasses (default sets surface receptors to zero)

    Random susceptibility is implemented by assigning a probability of a cell being initialized as unsusceptible
    according to a fractional value of unsusceptible cells.

    To track susceptibility, set the attribute `track_susc` to True.
    """

    unique_key = random_susceptibility_steppable_key

    def __init__(self, frequency=1):
        nCoVSteppableBase.__init__(self, frequency)

        self._track_susc = SusceptibilityModelInputs.track_susc

        self._epithelial_types = {}

        # Initialize default data
        self.append_epithelial_type(uninfected_type_name, SusceptibilityModelInputs.frac_not_susc)

    def start(self):
        """
        Called once when initializing a simulation
        """
        if self.track_susc:
            self.track_cell_level_scalar_attribute(field_name='Susceptible', attribute_name=vs_state_key)

    def make_unsusceptible(self, _cell) -> bool:
        """
        Try to make cell not susceptible to infection; can be overwritten in subclasses

        :param _cell: cell to try to make not susceptible
        :return: True if cell is made not susceptible; False otherwise
        """
        if unbound_receptors_cellg_key not in _cell.dict.keys():
            raise AttributeError('RandomSusceptibilitySteppable requires surface receptors')
        if _cell.dict[unbound_receptors_cellg_key] > 0:
            _cell.dict[unbound_receptors_cellg_key] = 0
            if self.track_susc:
                _cell.dict[vs_state_key] = False
            return True
        return False

    def append_epithelial_type(self, _name: str, _val: float = SusceptibilityModelInputs.frac_not_susc):
        """
        Append epithelial cell type data

        :param _name: Name of epithelial cell type
        :param _val: Fraction not susceptible
        :return: None
        """
        if not 0 <= _val <= 1:
            raise ValueError('Fractional values must be in [0, 1]')
        self._epithelial_types[_name] = _val

    def remove_epithelial_type(self, _name: str):
        """
        Remove epithelial cell type data

        :param _name: Name of epithelial cell type
        :return: None
        """
        self._epithelial_types.pop(_name)

    @property
    def epithelial_data(self) -> dict:
        """
        Copy of epithelial cell type data
        """
        return self._epithelial_types.copy()

    @property
    def epithelial_type_ids(self):
        """
        Ids of epithelial cell types according to a cc3d simulation
        """
        return [getattr(self, x.upper()) for x in self._epithelial_types.keys()]

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback. If a cell's type is registered, then it becomes unsusceptible with a probability
        according to its registered fraction
        """
        if _new_cell.type in self.epithelial_type_ids:
            if self.track_susc:
                _new_cell.dict[vs_state_key] = True
            if random() <= self._epithelial_types[self.get_type_name_by_cell(_new_cell)]:
                self.make_unsusceptible(_new_cell)

    @property
    def track_susc(self) -> bool:
        """
        Flag to track susceptibility per cell
        """
        return self._track_susc

    @track_susc.setter
    def track_susc(self, _val: bool):
        self._track_susc = _val
