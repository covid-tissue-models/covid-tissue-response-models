###############################################################################################################
# To cite this model please use the following:
#
# T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
# James M. Osborne, Ellen M. Quardokus, James A. Glazier,
# "A modular framework for multiscale multicellular spatial modeling of viral infection, immune response and drug
# therapy timing and efficacy in epithelial tissues",
# bioRxiv 2020.04.27.064139
###############################################################################################################

import os
import sys

import math

from cc3d.core.PySteppables import *
import numpy as np

rng = np.random  # alias for random number generators (rng)

# Import project libraries and classes
sys.path.append(os.path.dirname(__file__))
from ViralInfectionVTMSteppableBasePy import *
import ViralInfectionVTMLib
from ViralInfectionVTMModelInputs import *

# Import toolkit
sys.path.append(os.path.dirname(os.path.dirname(__file__)))
from nCoVToolkit import nCoVUtils


class CellsInitializerSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Initializes epithelial sheet and an initial immune cell population
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

    def start(self):
        self.get_xml_element('virus_dc').cdata = virus_dc
        self.get_xml_element('virus_decay').cdata = virus_decay

        for x in range(0, self.dim.x, int(cell_diameter)):
            for y in range(0, self.dim.y, int(cell_diameter)):
                cell = self.new_uninfected_cell_in_time()
                self.cellField[x:x + int(cell_diameter), y:y + int(cell_diameter), 0] = cell
                cell.dict[ViralInfectionVTMLib.vrl_key] = False
                ViralInfectionVTMLib.reset_viral_replication_variables(cell=cell)
                cell.dict['Receptors'] = initial_unbound_receptors
                self.load_viral_replication_model(cell=cell, vr_step_size=vr_step_size,
                                                  unpacking_rate=unpacking_rate,
                                                  replicating_rate=replicating_rate,
                                                  r_half=r_half,
                                                  translating_rate=translating_rate,
                                                  packing_rate=packing_rate)

        # Infect a cell
        cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
        cell.dict['Unpacking'] = 1.0
        cell.type = self.INFECTED

        self.load_viral_replication_model(cell=cell, vr_step_size=vr_step_size,
                                          unpacking_rate=unpacking_rate,
                                          replicating_rate=replicating_rate,
                                          r_half=r_half,
                                          translating_rate=translating_rate,
                                          packing_rate=packing_rate,
                                          secretion_rate=secretion_rate)

        cell.dict['ck_production'] = max_ck_secrete_infect

        for iteration in range(int(initial_immune_seeding)):
            cell = True
            while cell:
                xi = np.random.randint(0, self.dim.x - 2 * cell_diameter)
                yi = np.random.randint(0, self.dim.y - 2 * cell_diameter)
                for x in range(xi, xi + int(cell_diameter)):
                    for y in range(yi, yi + int(cell_diameter)):
                        cell = self.cell_field[x, y, 1]
                        break
                cell = False
            cell = self.new_immune_cell_in_time(ck_production=max_ck_secrete_im, ck_consumption=max_ck_consume)
            self.cell_field[x:x + int(cell_diameter), y:y + int(cell_diameter), 1] = cell
            cell.targetVolume = cell_volume
            cell.lambdaVolume = cell_volume
            cell.dict['activated'] = False  # flag for immune cell being naive or activated
            # cyttokine params
            cell.dict['ck_production'] = max_ck_secrete_im
            cell.dict['ck_consumption'] = max_ck_consume
            cell.dict['tot_ck_upt'] = 0


# TODO Add actual uptake from the field based on discussion with James
class ViralReplicationSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral replication module
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        self.plot_win = None

        # Reference to SimDataSteppable
        self.simdata_steppable = None

    def start(self):

        # Load model
        options = {'relative': 1e-10, 'absolute': 1e-12}
        self.set_sbml_global_options(options)

    def step(self, mcs):
        if self.simdata_steppable is None:
            self.simdata_steppable: SimDataSteppable = \
                self.shared_steppable_vars[ViralInfectionVTMLib.simdata_steppable_key]

        # Sample state of cell at center of domain (first infected cell)
        cell = self.cellField[self.dim.x / 2, self.dim.y / 2, 0]
        self.simdata_steppable.set_vrm_tracked_cell(cell=cell)

        # Do viral model
        for cell in self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING):
            # Step the model for this cell
            ViralInfectionVTMLib.step_sbml_model_cell(cell=cell)
            # Pack state variables into cell dictionary
            ViralInfectionVTMLib.pack_viral_replication_variables(cell=cell)

            # Test for infection secretion
            if cell.dict['Assembled'] > cell_infection_threshold:
                cell.type = self.INFECTEDSECRETING
                ViralInfectionVTMLib.enable_viral_secretion(cell=cell, secretion_rate=secretion_rate)

                # cytokine params
                cell.dict['ck_production'] = max_ck_secrete_infect

            # Test for cell death
            if cell.type == self.INFECTEDSECRETING and \
                    np.random.random() < nCoVUtils.hill_equation(cell.dict['Assembled'],
                                                                 diss_coeff_uptake_apo,
                                                                 hill_coeff_uptake_apo):
                self.kill_cell(cell=cell)


class ViralInternalizationSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral internalization module
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

    def start(self):
        # Post reference to self
        self.shared_steppable_vars[ViralInfectionVTMLib.vim_steppable_key] = self

    def step(self, mcs):
        pass

    def do_cell_internalization(self, cell, viral_amount_com):
        if cell.dict['Receptors'] == 0:
            return False, 0.0

        _k = kon * cell.volume / koff
        diss_coeff_uptake_pr = math.sqrt(initial_unbound_receptors / 2.0 / _k / cell.dict['Receptors'])
        uptake_probability = nCoVUtils.hill_equation(viral_amount_com,
                                                     diss_coeff_uptake_pr,
                                                     hill_coeff_uptake_pr)

        cell_does_uptake = np.random.rand() < uptake_probability
        uptake_amount = s_to_mcs / rate_coeff_uptake_pr * uptake_probability

        if cell_does_uptake and cell.type == self.UNINFECTED:
            cell.type = self.INFECTED
            cell.dict['ck_production'] = max_ck_secrete_infect
            self.load_viral_replication_model(cell=cell, vr_step_size=vr_step_size,
                                              unpacking_rate=unpacking_rate,
                                              replicating_rate=replicating_rate,
                                              r_half=r_half,
                                              translating_rate=translating_rate,
                                              packing_rate=packing_rate,
                                              secretion_rate=secretion_rate)

        return cell_does_uptake, uptake_amount

    def update_cell_receptors(self, cell, receptors_increment):
        cell.dict['Receptors'] = max(cell.dict['Receptors'] + receptors_increment, 0.0)


class ViralSecretionSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral release module
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        # Reference to ViralInternalizationSteppable
        self.vim_steppable = None

    def start(self):
        if track_model_variables:
            self.track_cell_level_scalar_attribute(field_name='Uptake', attribute_name='Uptake')
            self.track_cell_level_scalar_attribute(field_name='Assembled', attribute_name='Assembled')
            self.track_cell_level_scalar_attribute(field_name='Unpacking', attribute_name='Unpacking')
            self.track_cell_level_scalar_attribute(field_name='Replicating', attribute_name='Replicating')
            self.track_cell_level_scalar_attribute(field_name='Uptake', attribute_name='Uptake')
            self.track_cell_level_scalar_attribute(field_name='Secretion', attribute_name='Secretion')

    def step(self, mcs):
        if self.vim_steppable is None:
            self.vim_steppable: ViralInternalizationSteppable = \
                self.shared_steppable_vars[ViralInfectionVTMLib.vim_steppable_key]

        secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.INFECTEDSECRETING):

            # Evaluate probability of cell uptake of viral particles from environment
            # If cell isn't infected, it changes type to infected here if uptake occurs
            viral_amount_com = self.field.Virus[cell.xCOM, cell.yCOM, cell.zCOM] * cell.volume
            cell_does_uptake, uptake_amount = self.vim_steppable.do_cell_internalization(cell, viral_amount_com)
            if cell_does_uptake:
                uptake = secretor.uptakeInsideCellTotalCount(cell, 1E12, uptake_amount / cell.volume)
                cell.dict['Uptake'] = abs(uptake.tot_amount)
                self.vim_steppable.update_cell_receptors(cell=cell, receptors_increment=-cell.dict['Uptake'] * s_to_mcs)
                ViralInfectionVTMLib.set_viral_replication_cell_uptake(cell=cell, uptake=cell.dict['Uptake'])

            if cell.type == self.INFECTEDSECRETING:
                sec_amount = ViralInfectionVTMLib.get_viral_replication_cell_secretion(cell=cell)
                secretor.secreteInsideCellTotalCount(cell, sec_amount / cell.volume)


class ImmuneCellKillingSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell direct cytotoxicity and bystander effect module
    """

    def step(self, mcs):
        killed_cells = []
        for cell in self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.IMMUNECELL:
                        self.kill_cell(cell=cell)
                        killed_cells.append(cell)

        # Bystander Effect
        for cell in killed_cells:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type in [self.INFECTED, self.INFECTEDSECRETING, self.UNINFECTED]:
                        p_bystander_effect = np.random.random()
                        if p_bystander_effect < bystander_effect:
                            self.kill_cell(cell=neighbor)


class ChemotaxisSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell chemotaxis module
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list_by_type(self.IMMUNECELL):

            cd = self.chemotaxisPlugin.addChemotaxisData(cell, "cytokine")
            if cell.dict['activated']:
                cd.setLambda(lamda_chemotaxis)
            else:
                cd.setLambda(0.0)
            cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

    def step(self, mcs):
        field = self.field.Virus
        for cell in self.cell_list_by_type(self.IMMUNECELL):

            cd = self.chemotaxisPlugin.getChemotaxisData(cell, "cytokine")
            concentration = field[cell.xCOM, cell.yCOM, 1]
            constant = lamda_chemotaxis
            l = constant / (1.0 + concentration)
            if cell.dict['activated']:
                cd.setLambda(l)
            else:
                cd.setLambda(0)


class ImmuneCellSeedingSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell seeding and removal of immune cell recruitment module
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        # Reference to ImmuneResponseSteppable
        self.ir_steppable = None

    def step(self, mcs):
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = \
                self.shared_steppable_vars[ViralInfectionVTMLib.ir_steppable_key]

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            p_immune_dying = np.random.random()
            if p_immune_dying < self.ir_steppable.get_immune_removal_prob():
                cell.targetVolume = 0.0

        p_immune_seeding = np.random.random()
        if p_immune_seeding < self.ir_steppable.get_immune_seeding_prob():
            open_space = True
            viral_concentration = 0
            for iteration in range(10):
                radius = 10
                length = 0
                while length <= radius:
                    xi = np.random.randint(0, self.dim.x - 2 * cell_diameter)
                    yi = np.random.randint(0, self.dim.y - 2 * cell_diameter)
                    length = np.sqrt((self.dim.x // 2 - xi) ** 2 + (self.dim.y // 2 - yi) ** 2)
                for x in range(xi, xi + int(cell_diameter)):
                    for y in range(yi, yi + int(cell_diameter)):
                        cell = self.cellField[x, y, 1]
                        if cell:
                            open_space = False
                            break
                if open_space:
                    concentration_iteration = self.field.Virus[xi, yi, 1]
                    if concentration_iteration >= viral_concentration:
                        viral_concentration = concentration_iteration
                        x_seed = xi
                        y_seed = yi
            if open_space:
                cell = self.new_immune_cell_in_time(ck_production=max_ck_secrete_im, ck_consumption=max_ck_consume)

                self.cellField[x_seed:x_seed + int(cell_diameter), y_seed:y_seed + int(cell_diameter), 1] = cell
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "cytokine")
                if cell.dict['activated']:
                    cd.setLambda(lamda_chemotaxis)
                else:
                    cd.setLambda(0.0)

                cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
                cell.targetVolume = cell_volume
                cell.lambdaVolume = cell_volume


class SimDataSteppable(SteppableBasePy):
    """
    Plots/writes simulation data of interest
    """

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

        self.vrm_data_win = None
        self.vrm_data_path = None
        # The viral replication model of this cell is tracked and plotted/recorded
        self.vrm_tracked_cell = None

        self.vim_data_win = None
        self.vim_data_path = None

        self.pop_data_win = None
        self.pop_data_path = None

        self.med_diff_data_win = None
        self.med_diff_data_path = None

        self.ir_data_win = None
        self.ir_data_path = None

        self.spat_data_win = None
        self.spat_data_path = None

        self.plot_vrm_data = plot_vrm_data_freq > 0
        self.write_vrm_data = write_vrm_data_freq > 0

        self.plot_vim_data = plot_vim_data_freq > 0
        self.write_vim_data = write_vim_data_freq > 0

        self.plot_pop_data = plot_pop_data_freq > 0
        self.write_pop_data = write_pop_data_freq > 0

        self.plot_med_diff_data = plot_med_diff_data_freq > 0
        self.write_med_diff_data = write_med_diff_data_freq > 0
        self.med_diff_key = "MedDiff"

        self.plot_ir_data = plot_ir_data_freq > 0
        self.write_ir_data = write_ir_data_freq > 0
        self.ir_key = "ImmuneResp"
        self.ir_steppable = None

        self.plot_spat_data = plot_spat_data_freq > 0
        self.write_spat_data = write_spat_data_freq > 0

        # Origin of infection point; if more than one cell is first detected, then measure the mean COM
        # If first infection is far from center of domain, then measurements of infection front probably won't
        # be very useful
        self.init_infect_pt = None

    def start(self):
        # Post reference to self
        self.shared_steppable_vars[ViralInfectionVTMLib.simdata_steppable_key] = self

        if self.plot_vrm_data:
            self.vrm_data_win = self.add_new_plot_window(title='VRM',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Variables', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=False,
                                                         config_options={'legend': True})

            self.vrm_data_win.add_plot("U", style='Dots', color='blue', size=5)
            self.vrm_data_win.add_plot("R", style='Dots', color='orange', size=5)
            self.vrm_data_win.add_plot("P", style='Dots', color='green', size=5)
            self.vrm_data_win.add_plot("A", style='Dots', color='red', size=5)
            self.vrm_data_win.add_plot("Uptake", style='Dots', color='yellow', size=5)
            self.vrm_data_win.add_plot("Secretion", style='Dots', color='white', size=5)

        if self.plot_vim_data:
            self.vim_data_win = self.add_new_plot_window(title='VIM',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Variables', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=False,
                                                         config_options={'legend': True})

            self.vim_data_win.add_plot("R", style='Dots', color='orange', size=5)

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
            self.pop_data_win.add_plot("InfectedSecreting", style='Dots', color='green', size=5)
            self.pop_data_win.add_plot("Dying", style='Dots', color='yellow', size=5)
            self.pop_data_win.add_plot("ImmuneCell", style='Dots', color='white', size=5)
            self.pop_data_win.add_plot("ImmuneCellActivated", style='Dots', color='purple', size=5)

        if self.plot_med_diff_data:
            self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                              x_axis_title='MCS',
                                                              y_axis_title='Number of diffusive species per volume',
                                                              x_scale_type='linear',
                                                              y_scale_type='log',
                                                              grid=True,
                                                              config_options={'legend': True})

            self.med_diff_data_win.add_plot("MedViral", style='Dots', color='red', size=5)
            self.med_diff_data_win.add_plot("MedCyt", style='Dots', color='blue', size=5)
            self.med_diff_data_win.add_plot("MedOxi", style='Dots', color='green', size=5)

        if self.plot_ir_data:
            self.ir_data_win = self.add_new_plot_window(title='Immune Response Model',
                                                        x_axis_title='MCS',
                                                        y_axis_title='State variable S',
                                                        x_scale_type='linear',
                                                        y_scale_type='linear',
                                                        grid=True)

            self.ir_data_win.add_plot(self.ir_key, style='Dots', color='red', size=5)

        if self.plot_spat_data:
            self.spat_data_win = self.add_new_plot_window(title='Spatial data',
                                                          x_axis_title='MCS',
                                                          y_axis_title='',
                                                          x_scale_type='linear',
                                                          y_scale_type='linear',
                                                          grid=True,
                                                          config_options={'legend': True})

            self.spat_data_win.add_plot("DeathComp", style='Dots', color='red', size=5)
            self.spat_data_win.add_plot("InfectDist", style='Dots', color='blue', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_vrm_data:
                self.vrm_data_path = Path(self.output_dir).joinpath('vrm_data.dat')
                with open(self.vrm_data_path, 'w'):
                    pass

            if self.write_vim_data:
                self.vim_data_path = Path(self.output_dir).joinpath('vim_data.dat')
                with open(self.vim_data_path, 'w'):
                    pass

            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath('pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

            if self.write_med_diff_data:
                self.med_diff_data_path = Path(self.output_dir).joinpath('med_diff_data.dat')
                with open(self.med_diff_data_path, 'w'):
                    pass

            if self.write_ir_data:
                self.ir_data_path = Path(self.output_dir).joinpath('ir_data.dat')
                with open(self.ir_data_path, 'w'):
                    pass

            if self.write_spat_data:
                self.spat_data_path = Path(self.output_dir).joinpath('spat_data.dat')
                with open(self.spat_data_path, 'w'):
                    pass

    def step(self, mcs):

        plot_pop_data = self.plot_pop_data and mcs % plot_pop_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % plot_med_diff_data_freq == 0
        plot_ir_data = self.plot_ir_data and mcs % plot_ir_data_freq == 0
        plot_vrm_data = self.plot_vrm_data and mcs % plot_vrm_data_freq == 0
        plot_vim_data = self.plot_vim_data and mcs % plot_vim_data_freq == 0
        plot_spat_data = self.plot_spat_data and mcs % plot_spat_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % write_pop_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % write_med_diff_data_freq == 0
            write_ir_data = self.write_ir_data and mcs % write_ir_data_freq == 0
            write_vrm_data = self.write_vrm_data and mcs % write_vrm_data_freq == 0
            write_vim_data = self.write_vim_data and mcs % write_vim_data_freq == 0
            write_spat_data = self.write_spat_data and mcs % write_spat_data_freq == 0
        else:
            write_pop_data = False
            write_med_diff_data = False
            write_ir_data = False
            write_vrm_data = False
            write_vim_data = False
            write_spat_data = False

        if self.vrm_tracked_cell is not None and (plot_vrm_data or write_vrm_data):
            if plot_vrm_data:
                self.vrm_data_win.add_data_point("U", mcs, self.vrm_tracked_cell.dict['Unpacking'])
                self.vrm_data_win.add_data_point("R", mcs, self.vrm_tracked_cell.dict['Replicating'])
                self.vrm_data_win.add_data_point("P", mcs, self.vrm_tracked_cell.dict['Packing'])
                self.vrm_data_win.add_data_point("A", mcs, self.vrm_tracked_cell.dict['Assembled'])
                self.vrm_data_win.add_data_point("Uptake", mcs, self.vrm_tracked_cell.dict['Uptake'])
                self.vrm_data_win.add_data_point("Secretion", mcs, self.vrm_tracked_cell.dict['Secretion'])

            if write_vrm_data:
                with open(self.vrm_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}, {}, {}, {}, {}\n'.format(mcs,
                                                                         self.vrm_tracked_cell.id,
                                                                         self.vrm_tracked_cell.dict['Unpacking'],
                                                                         self.vrm_tracked_cell.dict['Replicating'],
                                                                         self.vrm_tracked_cell.dict['Packing'],
                                                                         self.vrm_tracked_cell.dict['Assembled'],
                                                                         self.vrm_tracked_cell.dict['Uptake'],
                                                                         self.vrm_tracked_cell.dict['Secretion']))

        if self.vrm_tracked_cell is not None and (plot_vim_data or write_vim_data):
            if plot_vim_data:
                self.vim_data_win.add_data_point("R", mcs, self.vrm_tracked_cell.dict['Receptors'])

            if write_vim_data:
                with open(self.vim_data_path, 'a') as fout:
                    fout.write('{}, {}, {}\n'.format(mcs,
                                                     self.vrm_tracked_cell.id,
                                                     self.vrm_tracked_cell.dict['Receptors']))

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.UNINFECTED))
            num_cells_infected = len(self.cell_list_by_type(self.INFECTED))
            num_cells_infectedsecreting = len(self.cell_list_by_type(self.INFECTEDSECRETING))
            num_cells_dying = len(self.cell_list_by_type(self.DYING))
            num_cells_immune = len(self.cell_list_by_type(self.IMMUNECELL))
            num_cells_immune_act = len([c for c in self.cell_list_by_type(self.IMMUNECELL) if c.dict['activated']])

            # Plot population data plot if requested
            if plot_pop_data:
                if num_cells_uninfected > 0:
                    self.pop_data_win.add_data_point('Uninfected', mcs, num_cells_uninfected)
                if num_cells_infected > 0:
                    self.pop_data_win.add_data_point('Infected', mcs, num_cells_infected)
                if num_cells_infectedsecreting > 0:
                    self.pop_data_win.add_data_point('InfectedSecreting', mcs, num_cells_infectedsecreting)
                if num_cells_dying > 0:
                    self.pop_data_win.add_data_point('Dying', mcs, num_cells_dying)
                if num_cells_immune > 0:
                    self.pop_data_win.add_data_point('ImmuneCell', mcs, num_cells_immune)
                if num_cells_immune_act > 0:
                    self.pop_data_win.add_data_point('ImmuneCellActivated', mcs, num_cells_immune_act)

            # Write population data to file if requested
            if write_pop_data:
                with open(self.pop_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}, {}, {}, {}\n'.format(mcs,
                                                                     num_cells_uninfected,
                                                                     num_cells_infected,
                                                                     num_cells_infectedsecreting,
                                                                     num_cells_dying,
                                                                     num_cells_immune,
                                                                     num_cells_immune_act))

        if plot_med_diff_data or write_med_diff_data:

            # Gather total diffusive amounts
            med_viral_total = 0.0
            med_cyt_total = 0.0
            med_oxi_total = 0.0
            for x, y, z in self.every_pixel():
                med_viral_total += self.field.Virus[x, y, z]
                med_cyt_total += self.field.cytokine[x, y, z]
                med_oxi_total += self.field.oxidator[x, y, z]

            # Plot total diffusive viral amount if requested
            if plot_med_diff_data:
                if med_viral_total > 0:
                    self.med_diff_data_win.add_data_point("MedViral", mcs, med_viral_total)
                if med_cyt_total > 0:
                    self.med_diff_data_win.add_data_point("MedCyt", mcs, med_cyt_total)
                if med_oxi_total > 0:
                    self.med_diff_data_win.add_data_point("MedOxi", mcs, med_oxi_total)

            # Write total diffusive viral amount if requested
            if write_med_diff_data:
                with open(self.med_diff_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}\n'.format(mcs, med_viral_total, med_cyt_total, med_oxi_total))

        if plot_ir_data or write_ir_data:
            if self.ir_steppable is None:
                if self.ir_steppable is None:
                    self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[
                        ViralInfectionVTMLib.ir_steppable_key]

            s_val = self.ir_steppable.get_state_variable_val()

            # Plot state variable S if requested
            if plot_ir_data:
                self.ir_data_win.add_data_point(self.ir_key, mcs, s_val)

            # Write state variable S if requested
            if write_ir_data:
                with open(self.ir_data_path, 'a') as fout:
                    fout.write('{}, {}\n'.format(mcs, s_val))

        if plot_spat_data or write_spat_data:
            # Calculate compactness of dead cell area as total surface area of intefaces between dying and non-dying
            # types in epithelial sheet divided by total volume of dying types
            dead_srf = 0
            dead_vol = 0
            dying_cell_list = self.cell_list_by_type(self.DYING)
            if not dying_cell_list:
                dead_comp = 0
            else:
                for cell in dying_cell_list:
                    dead_vol += cell.volume
                    for neighbor, common_srf in self.get_cell_neighbor_data_list(cell):
                        if neighbor is not None and neighbor.type in [self.UNINFECTED,
                                                                      self.INFECTED,
                                                                      self.INFECTEDSECRETING]:
                            dead_srf += common_srf

                dead_comp = dead_srf / dead_vol

            # Calculate infection front: max. distance from initial point of infection to all infected cells
            # If no infected cells, distance is -1
            max_infect_dist = -1
            if self.init_infect_pt is None:
                infected_cell_list = self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING)
                num_cells_infected = len(infected_cell_list)
                if num_cells_infected > 0:
                    self.init_infect_pt = [0, 0, 0]
                    for cell in infected_cell_list:
                        self.init_infect_pt[0] += cell.xCOM
                        self.init_infect_pt[1] += cell.yCOM

                    self.init_infect_pt[0] /= num_cells_infected
                    self.init_infect_pt[1] /= num_cells_infected

            if self.init_infect_pt is not None:
                for cell in self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING):
                    dx = cell.xCOM - self.init_infect_pt[0]
                    dy = cell.yCOM - self.init_infect_pt[1]
                    max_infect_dist = max(max_infect_dist, math.sqrt(dx * dx + dy * dy))

            # Plot spatial data if requested
            #   Infection distance is normalized by average lattice dimension
            if plot_spat_data:
                self.spat_data_win.add_data_point("DeathComp", mcs, dead_comp)
                if max_infect_dist > 0:
                    max_infect_dist_norm = max_infect_dist / ((self.dim.x + self.dim.y) / 2.0)
                    self.spat_data_win.add_data_point("InfectDist", mcs, max_infect_dist_norm)

            # Write spatial data if requested
            if write_spat_data:
                with open(self.spat_data_path, 'a') as fout:
                    fout.write('{}, {}, {}\n'.format(mcs, dead_comp, max_infect_dist))

    def set_vrm_tracked_cell(self, cell):
        self.vrm_tracked_cell = cell


class CytokineProductionAbsorptionSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements cytokine production/secretion and immune cell activation module
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        if track_model_variables:
            self.track_cell_level_scalar_attribute(field_name='activated', attribute_name='activated')
        self.ck_secretor = None
        self.virus_secretor = None

        # Reference to ImmuneResponseSteppable
        self.ir_steppable = None

    def start(self):
        # cytokine diff parameters
        self.get_xml_element('cytokine_dc').cdata = cytokine_dc
        self.get_xml_element('cytokine_decay').cdata = cytokine_field_decay  # no "natural" decay, only consumption
        # and "leakage" outside of simulation laticce

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            # cytokine production/uptake parameters for immune cells

            cell.dict['ck_production'] = max_ck_secrete_im
            cell.dict['ck_consumption'] = max_ck_consume

        for cell in self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING):
            cell.dict['ck_production'] = max_ck_secrete_infect

        self.ck_secretor = self.get_field_secretor("cytokine")
        self.virus_secretor = self.get_field_secretor("Virus")

    def step(self, mcs):
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = \
                self.shared_steppable_vars[ViralInfectionVTMLib.ir_steppable_key]

        # Track the total amount added and subtracted to the cytokine field
        total_ck_inc = 0.0

        for cell in self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING):
            viral_load = ViralInfectionVTMLib.get_assembled_viral_load_inside_cell(cell, vr_step_size)
            produced = cell.dict['ck_production'] * nCoVUtils.hill_equation(viral_load, ec50_infecte_ck_prod, 2)
            res = self.ck_secretor.secreteInsideCellTotalCount(cell, produced / cell.volume)
            total_ck_inc += res.tot_amount

        for cell in self.cell_list_by_type(self.IMMUNECELL):

            self.virus_secretor.uptakeInsideCellTotalCount(cell, cell.dict['ck_consumption'] / cell.volume, 0.1)

            up_res = self.ck_secretor.uptakeInsideCellTotalCount(cell,
                                                                 cell.dict['ck_consumption'] / cell.volume, 0.1)
            # decay seen ck
            cell.dict['tot_ck_upt'] *= ck_memory_immune

            # uptake ck

            cell.dict['tot_ck_upt'] -= up_res.tot_amount  # from POV of secretion uptake is negative
            total_ck_inc += up_res.tot_amount
            p_activate = nCoVUtils.hill_equation(cell.dict['tot_ck_upt'], EC50_ck_immune, 2)

            if rng.uniform() < p_activate and not cell.dict['activated']:

                cell.dict['activated'] = True
                cell.dict['time_activation'] = mcs
            elif (cell.dict['activated']
                  and mcs - cell.dict['time_activation'] > minimum_activated_time):
                cell.dict['activated'] = False
                cell.dict['time_activation'] = - 99

            if cell.dict['activated']:
                seen_field = self.total_seen_field(self.field.cytokine, cell)
                produced = cell.dict['ck_production'] * nCoVUtils.hill_equation(seen_field, 100, 1)
                sec_res = self.ck_secretor.secreteInsideCellTotalCount(cell, produced / cell.volume)

                total_ck_inc += sec_res.tot_amount

        self.ir_steppable.increment_total_cytokine_count(total_ck_inc)


class ImmuneRecruitmentSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell recruitment module
    Note that total cytokine is currently tracked elsewhere by counting uptake and secretion, and by applying the
    field decay rate applied to the cytokine field. This is only relevant for periodic and zero-flux boundary conditions
    for the diffusive cytokine field.
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        # Reference to solver
        self.__rr = None

        # Running value of total cytokine; to be updated externally through accessor
        self.__total_cytokine = 0.0

        self.__ck_decay = 0.0

    def start(self):
        self.__ck_decay = float(self.get_xml_element('cytokine_decay').cdata)

        # Post reference to self
        self.shared_steppable_vars[ViralInfectionVTMLib.ir_steppable_key] = self

        # Initialize model
        self.__init_fresh_recruitment_model()

    def step(self, mcs):

        # Update total count of immune cells
        num_immune_cells = len(self.cell_list_by_type(self.IMMUNECELL))

        # Apply consumption / transmission decay to running total
        total_cytokine_decayed = self.__total_cytokine * self.__ck_decay
        self.__total_cytokine -= total_cytokine_decayed

        # Update model
        total_cytokine_transmitted = ir_transmission_coeff * total_cytokine_decayed
        self.update_running_recruitment_model(num_immune_cells, total_cytokine_transmitted)

    def finish(self):
        pass

    def __init_fresh_recruitment_model(self):
        # Generate solver instance
        model_string = ViralInfectionVTMLib.immune_recruitment_model_string(ir_add_coeff,
                                                                            ir_subtract_coeff,
                                                                            ir_delay_coeff,
                                                                            ir_decay_coeff)
        self.add_free_floating_antimony(model_string=model_string,
                                        model_name=ViralInfectionVTMLib.ir_model_name,
                                        step_size=vr_step_size)

        # Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == ViralInfectionVTMLib.ir_model_name:
                self.__rr = rr

    def update_running_recruitment_model(self, num_immune_cells, total_cytokine):
        self.__rr['numImmuneCells'] = num_immune_cells
        self.__rr['totalCytokine'] = total_cytokine
        self.__rr.timestep()

    def get_state_variable_val(self):
        return self.__rr['S']

    def get_immune_seeding_prob(self):
        """
        Returns probability of immune cell seeding due to local and global recruitment
        Probability is only non-zero if the state variable *S* is positive, in which case
        the probability is an error function of *S*
        :return: probability of immune cell seeding due to local and global recruitment
        """
        s_val = self.get_state_variable_val()
        if s_val < 0:
            return 0.0
        else:
            return math.erf(ir_prob_scaling_factor * s_val)

    def get_immune_removal_prob(self):
        """
        Returns probability of immune cell removal due to local and global recruitment
        Probability is only non-zero if the state variable *S* is negative, in which case
        the probability is an error function of *S*
        :return: probability of immune cell removal due to local and global recruitment
        """
        s_val = self.get_state_variable_val()
        if s_val > 0:
            return 0.0
        else:
            return math.erf(- ir_prob_scaling_factor * s_val)

    def increment_total_cytokine_count(self, _inc_amount):
        """
        Method to efficient maintain total count of cytokine
        :param _inc_amount: total amount increment
        :return: None
        """
        self.__total_cytokine = max(0.0, self.__total_cytokine + _inc_amount)


class oxidationAgentModelSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell oxidizing agent cytotoxicity module
    """
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)
        if track_model_variables:
            self.track_cell_level_scalar_attribute(field_name='oxi_killed', attribute_name='oxi_killed')
        self.oxi_secretor = None

    def start(self):
        self.get_xml_element('oxi_dc').cdata = oxi_dc
        self.get_xml_element('oxi_decay').cdata = oxi_decay

        self.oxi_secretor = self.get_field_secretor("oxidator")

    def step(self, mcs):

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            if cell.dict['activated']:
                seen_field = self.total_seen_field(self.field.cytokine, cell)
                if seen_field > oxi_sec_thr:
                    oxi_sec = self.oxi_secretor.secreteInsideCellTotalCount(cell, max_oxi_secrete / cell.volume)

        for cell in self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.INFECTEDSECRETING):

            seen_field = self.total_seen_field(self.field.oxidator, cell)
            if seen_field >= oxi_death_thr:
                self.kill_cell(cell=cell)
                cell.dict['oxi_killed'] = True

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return
