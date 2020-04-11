# This is a library of steppable classes for the coronavirus viral infection modeling project using CompuCell3D
# by the Biocomplexity Institute at Indiana University using CompuCell3D

import os
import sys

from cc3d.cpp import CompuCell
import numpy as np

# Set this to True for local references when developing; False when running
__dev_mode__ = False

# This is just a gentle reminder to turn off debug mode
assert not __dev_mode__, "Don't forget to set to run!"

if __dev_mode__:
    # Import project libraries
    from . import CoronavirusLib

    # Import toolkit
    from ..nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
    from ..nCoVToolkit import nCoVUtils

else:
    # Import project libraries
    sys.path.append(os.path.dirname(__file__))
    import CoronavirusLib

    # Import toolkit
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
    from nCoVToolkit import nCoVUtils


class CoronavirusSteppableBasePy(nCoVSteppableBase):

    def __init__(self, frequency=1):
        nCoVSteppableBase.__init__(self, frequency)

    def start(self):
        pass

    def step(self, mcs):
        pass

    def finish(self):
        pass

    def load_viral_replication_model(self, cell, vr_step_size, unpacking_rate=0, replicating_rate=0,
                                     translating_rate=0, packing_rate=0, secretion_rate=0):
        """
        Loads viral replication model for a cell; initial values of state model are extract from cell.dict
        :param cell: cell for which the viral replication model is loaded
        :param vr_step_size: Antimony/SBML model step size
        :param unpacking_rate: model unpacking rate
        :param replicating_rate: model replicating rate
        :param translating_rate: model translating rate
        :param packing_rate: model packing rate
        :param secretion_rate: model secretion rate
        :return: None
        """
        if cell.dict[CoronavirusLib.vrl_key]:
            self.delete_sbml_from_cell(CoronavirusLib.vr_model_name, cell)

        model_string = CoronavirusLib.viral_replication_model_string(
            unpacking_rate, replicating_rate, translating_rate, packing_rate, secretion_rate,
            cell.dict['Unpacking'], cell.dict['Replicating'], cell.dict['Packing'], cell.dict['Assembled'],
            cell.dict['Uptake'])
        self.add_antimony_to_cell(model_string=model_string,
                                  model_name=CoronavirusLib.vr_model_name,
                                  cell=cell,
                                  step_size=vr_step_size)
        cell.dict[CoronavirusLib.vrl_key] = True
        CoronavirusLib.enable_viral_secretion(cell, cell.type == self.INFECTEDSECRETING)

    def cell_uptakes_virus(self, viral_field, cell, diss_coeff_uptake_pr, hill_coeff_uptake_pr, go_fast=True):
        """
        Calculates the probability of viral uptake from the environment as a function of local viral particle amount
        Returns true if cell uptakes virus
        Model development note: ACE2, TMPRSS2 effects may be well-implemented here in future work
        :param viral_field: environmental viral field
        :param cell: cell
        :param diss_coeff_uptake_pr: dissociation coefficient of Hill equation for probability function
        :param hill_coeff_uptake_pr: Hill coefficient of Hill equation for probability function
        :param go_fast: when True, a fast homogenized measurement is made; when False, it's slower but more accurate
        :return: True if cell uptakes; False if not
        """

        # Calculate total viral amount in cell's domain

        if go_fast:
            # Fast measurement
            cell_env_viral_val = viral_field[cell.xCOM, cell.yCOM, cell.zCOM] * cell.volume
        else:
            # Accurate measurement
            cell_env_viral_val = 0.0
            for ptd in self.get_cell_pixel_list(cell):
                cell_env_viral_val += viral_field[ptd.pixel.x, ptd.pixel.y, ptd.pixel.z]

        # Evaluate probability of uptake

        if cell_env_viral_val != 0:
            max_uptake_pr = nCoVUtils.hill_equation(val=cell_env_viral_val,
                                                    diss_cf=diss_coeff_uptake_pr,
                                                    hill_cf=hill_coeff_uptake_pr)
            return np.random.random() < max_uptake_pr
        else:
            return False

    def death_signal_binding(self, death_field, cell, diss_coeff_bind_pr, hill_coeff_bind_pr, trail):
        """
        Calculates the probability of extracellular death ligands binding to a cell receptor
        Returns true if ligand binds to receptor
        :param death_field: environmental viral field
        :param cell: cell
        :param diss_coeff_bind_pr: dissociation coefficient of Hill equation for probability function
        :param hill_coeff_bind_pr: Hill coefficient of Hill equation for probability function
        :param trail: percent of pre-assembled particles that leak out during packing
        :return: True if ligand binds to cell; False if not
        """
        cell_env_death_ligands = death_field[cell.xCOM, cell.yCOM, cell.zCOM] * cell.volume
        # print('cell_env_death_ligands=', cell_env_death_ligands)
        if cell_env_death_ligands != 0:
            max_death_ligands = nCoVUtils.hill_equation(val=cell_env_death_ligands,
                                                        diss_cf=diss_coeff_bind_pr,
                                                        hill_cf=hill_coeff_bind_pr)
            max_death_ligands = cell_env_death_ligands
            print(' max_death_ligands=',  max_death_ligands)
            result = np.random.random() < max_death_ligands * trail
            # if fraction of cells dying is too large, reduce trail to compensate
            # print('result=', result)
            return result
        else:
            return False

    def intrinsic_pathway(self, cell):
        """
        Calculates the probability that a cell will die via the intrinsic pathway as a rate based function of viral load
        :param cell: cell
        :return: True if probability of cell dying is met, False if not
        """
        # P = get_assembled_viral_load_inside_cell(cell=cell)



    def kill_cell(self, cell):
        """
        Model-specific cell death routines
        :param cell: cell to kill
        :return: None
        """
        cell.type = self.DYING
        CoronavirusLib.reset_viral_replication_variables(cell=cell)
        # Remove state model: no model for dead cell type
        self.remove_viral_replication_model(cell=cell)

    def remove_viral_replication_model(self, cell):
        """
        Removes viral replication model for a cell
        :param cell: cell for which to remove the viral replication model
        :return: None
        """
        if cell.dict[CoronavirusLib.vrl_key]:
            self.delete_sbml_from_cell(CoronavirusLib.vr_model_name, cell)
            cell.dict[CoronavirusLib.vrl_key] = False
