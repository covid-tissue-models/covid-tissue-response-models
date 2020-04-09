###############################################################################################################
# To cite this model please use the following:
#
# Josua Aponte-Serrano, T.J. Sego, Juliano F. Gianlupi, James A. Glazier,
# "Model of Viral Tissue Infection"
# https://github.com/covid-tissue-models/covid-tissue-response-models/tree/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM
###############################################################################################################

import os
import sys

from cc3d.core.PySteppables import *
from cc3d.cpp import CompuCell
import numpy as np

rng = np.random #alias for random number generators (rng)

# Set this to True for local references when developing; False when running
__dev_mode__ = False

# This is just a gentle reminder to turn off debug mode
assert not __dev_mode__, "Don't forget to set to run!"

if __dev_mode__:

    # Import project libraries and classes
    from .CoronavirusSteppableBasePy import *
    from . import CoronavirusLib

    # Import toolkit
    from ..nCoVToolkit import nCoVUtils

else:
    # Import project libraries and classes
    sys.path.append(os.path.dirname(__file__))
    from CoronavirusSteppableBasePy import *
    import CoronavirusLib

    # Import toolkit
    sys.path.append(os.path.dirname(os.path.dirname(__file__)))
    from nCoVToolkit import nCoVUtils

# Data control options
plot_pop_data_freq = 0  # Plot population data frequency (disable with 0)
write_pop_data_freq = 0  # Write population data to simulation directory frequency (disable with 0)
plot_med_viral_data_freq = 0  # Plot total diffusive viral amount frequency (disable with 0)
write_med_viral_data_freq = 0  # Write total diffusive viral amount frequency (disable with 0)

# Conversion Factors
s_to_mcs = 120.0  # s/mcs
um_to_lat_width = 4.0  # um/lattice_length

pmol_to_cc3d_au = 1e15

# Experimental Parameters
exp_cell_diameter = 12.0  # um

exp_replicating_rate = 1.0 / 20.0 * 1.0 / 60.0  # 1.0/20.0min * 1.0min/60.0s = 1.0/1200.0s
exp_translating_rate = exp_replicating_rate * 2.0
exp_unpacking_rate = exp_replicating_rate * 20.0
exp_packing_rate = exp_replicating_rate * 4.0
exp_secretion_rate = exp_replicating_rate * 4.0

exp_virus_dc = 10.0 / 100.0  # um^2/s

# cytokines:
# data from https://www.sciencedirect.com/science/article/pii/S1074761317300924 supplemental materials (A)
# and
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3433682/ (B)
cytoplasm_density = 6  # unit = [density of water](B)

exp_cytokine_dc_w = 100  # um^2/s; diffusion constant in water (A,B)
# the division by 10 is due to the small lattice size, as fiscussed with josh. we all should discuss this matter
exp_cytokine_dc_cyto = 16 / 10  # um^2/s; estimated diffusion constant in cytoplasm (B)
# ^ the  /10 is not experimental; added because of relative small area and because virus D is (or was) slowed down

exp_max_ck_diff_len = 100 # um  from (A); this diffusion length is due to uptake, I'm using it as a base
# for the decay due to ck leakage outside of the simulatted lattice.
# dl = sqrt(D/g) -> g = D/dl**2
exp_min_ck_decay = exp_cytokine_dc_cyto/(1.1*exp_max_ck_diff_len)**2

exp_max_cytokine_consumption = 1  # molecule / (cell second); maximum consumption of cytokine; actually a range [0.3,1] molecule / (cell second) (A)
exp_max_cytokine_immune_secretion = 10  # molecule / (cell second) (B)

exp_max_cytokine_consumption_mol = 3.5e-4  # pM/s
exp_max_cytokine_immune_secretion_mol = 3.5e-3  # pM/s

exp_EC50_cytokine_immune = 1  # pM from (B), it's a range from [1,50]pM

## tbd: try to find experimental data
minimum_activated_time_seconds = 60 * 60 # min * s/min


# =============================
# CompuCell Parameters
# cell
cell_diameter = exp_cell_diameter * 1.0 / um_to_lat_width
cell_volume = cell_diameter ** 2
# virus diffusion
virus_dc = exp_virus_dc * s_to_mcs / (um_to_lat_width ** 2)  # virus diffusion constant
virus_dl = cell_diameter * 3.0  # virus diffusion length
virus_decay = virus_dc / (virus_dl ** 2)  # virus decay rate
# virus intra-cellular
unpacking_rate = exp_unpacking_rate * s_to_mcs
replicating_rate = exp_replicating_rate * s_to_mcs
translating_rate = exp_translating_rate * s_to_mcs
packing_rate = exp_packing_rate * s_to_mcs
secretion_rate = exp_secretion_rate * s_to_mcs

# cytokine / immune activation

cytokine_dc = exp_cytokine_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)  # CK diff cst
# [1/g] = [s] -> [1/g]/s_to_mcs = [s]/[s/mcs] = [mcs]
# -> [g] = [1/s] -> [g]*[s_to_mcs] = [1/s]*[s/mcs] = [1/mcs]
cytokine_field_decay = exp_min_ck_decay * s_to_mcs
# pM = pmol/L = pmol/(10^15 um^3) = 10^-15 pmol/(um^3) = 10^-15 * um_to_lat_width^3 pmol/pixel
# pM/s = pM * s_to_mcs / MCS
max_ck_consume = exp_max_cytokine_consumption_mol * um_to_lat_width ** 3 * s_to_mcs * 1e-15 * pmol_to_cc3d_au  # cc3d_au/(pixel seconds)
max_ck_secrete_im = exp_max_cytokine_immune_secretion_mol * um_to_lat_width ** 3 * s_to_mcs * 1e-15 * pmol_to_cc3d_au  # * cc3d_au/(pixel seconds)
EC50_ck_immune = exp_EC50_cytokine_immune * um_to_lat_width ** 3 * 1e-15 * pmol_to_cc3d_au  # * cc3d_au/pixel
# ck_equilibrium = 1.5*EC50_ck_immune # equilibrium amount of ck in immune surface
ck_equilibrium = 2.1*EC50_ck_immune # equilibrium amount of ck in immune surface
ck_memory_immune = 1 - max_ck_consume/ck_equilibrium # decay therm for "seen" ck by immune

max_ck_secrete_infect = 10*max_ck_secrete_im
minimum_activated_time = minimum_activated_time_seconds/s_to_mcs # mcs

# Threshold at which cell infection is evaluated
cell_infection_threshold = 1.0
# Threshold at which cell death is evaluated
cell_death_threshold = 1.2
# Probability of survival of infected cell once cell_death_threshold is reached
survival_probability = 0.95

# Hill equation coefficients for probability of viral particle uptake from the environment
# Measurements are taken w.r.t. the total amount of viral particles in a cell's simulation subdomain
# dissociationt constant
diss_coeff_uptake_pr = 1.0
# Hill coefficient
hill_coeff_uptake_pr = 3.0

# Efficiency of viral uptake
relative_viral_uptake = 0.1

# Number of immune cells to seed at the beginning of the simulation
initial_immune_seeding = 0.0
# Rate for seeding of immune cells (constant)
immune_seeding_rate = 1.0 / 10.0
# Max dying rate of immune cells (actual rate is proportional to fraction of infected cells)
immunecell_dying_rate = 0.0 / 500.0
# Bystander effect
bystander_effect = 2.0/4.0
# Lambda Chemotaxis
# lamda_chemotaxis = 100.0
lamda_chemotaxis = 100.0/100.0

# Antimony/SBML model step size
vr_step_size = 1.0


class CellsInitializerSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
    """
    def __init__(self, frequency=1):
        CoronavirusSteppableBasePy.__init__(self, frequency)

    def start(self):
        self.get_xml_element('virus_dc').cdata = virus_dc
        self.get_xml_element('virus_decay').cdata = virus_decay

        for x in range(0, self.dim.x, int(cell_diameter)):
            for y in range(0, self.dim.y, int(cell_diameter)):
                cell = self.new_cell(self.UNINFECTED)
                self.cellField[x:x + int(cell_diameter), y:y + int(cell_diameter), 0] = cell
                cell.dict[CoronavirusLib.vrl_key] = False
                CoronavirusLib.reset_viral_replication_variables(cell=cell)
                cell.dict['Survived'] = False

        # Infect a cell
        cell = self.cell_field[self.dim.x // 2, self.dim.y // 2, 0]
        cell.dict['Unpacking'] = 1.0
        cell.type = self.INFECTED

        self.load_viral_replication_model(cell=cell, vr_step_size=vr_step_size,
                                          unpacking_rate=unpacking_rate,
                                          replicating_rate=replicating_rate,
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
                        cell = self.cellField[x, y, 1]
                        break
                cell = False
            cell = self.new_cell(self.IMMUNECELL)
            self.cellField[x:x + int(cell_diameter), y:y + int(cell_diameter), 1] = cell
            cell.targetVolume = cell_volume
            cell.lambdaVolume = cell_volume
            cell.dict['activated'] = False  # flag for immune cell being naive or activated
            # cyttokine params
            cell.dict['ck_production'] = max_ck_secrete_im  # TODO: replace secretion by hill
            cell.dict['ck_consumption'] = max_ck_consume  # TODO: replace by hill
            cell.dict['tot_ck_upt'] = 0


class Viral_ReplicationSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
    """
    def __init__(self, frequency=1):
        CoronavirusSteppableBasePy.__init__(self, frequency)

        self.plot_win = None

    def start(self):
        self.plot_win = self.add_new_plot_window(title='VRM',
                                                 x_axis_title='MonteCarlo Step (MCS)',
                                                 y_axis_title='Variables', x_scale_type='linear', y_scale_type='linear',
                                                 grid=False,
                                                 config_options={'legend': True})

        self.plot_win.add_plot("U", style='Dots', color='blue', size=5)
        self.plot_win.add_plot("R", style='Dots', color='orange', size=5)
        self.plot_win.add_plot("P", style='Dots', color='green', size=5)
        self.plot_win.add_plot("A", style='Dots', color='red', size=5)
        self.plot_win.add_plot("Uptake", style='Dots', color='yellow', size=5)
        self.plot_win.add_plot("Secretion", style='Dots', color='white', size=5)

        # Load model
        options = {'relative': 1e-10, 'absolute': 1e-12}
        self.set_sbml_global_options(options)

    def step(self, mcs):

        # Report rates to console
        print("Unpacking Rate = " + str(unpacking_rate))
        print("Replicating Rate = " + str(replicating_rate))
        print("Translating Rate = " + str(translating_rate))
        print("Packing Rate = " + str(packing_rate))
        print("Secretion Rate = " + str(secretion_rate))

        # Sample state of cell at center of domain (first infected cell)
        cell = self.cellField[self.dim.x / 2, self.dim.y / 2, 0]
        # Or sample state of cell near the first infected cell
        # cell = self.cellField[40, 45, 0]
        self.plot_win.add_data_point("U", mcs, cell.dict['Unpacking'])
        self.plot_win.add_data_point("R", mcs, cell.dict['Replicating'])
        self.plot_win.add_data_point("P", mcs, cell.dict['Packing'])
        self.plot_win.add_data_point("A", mcs, cell.dict['Assembled'])
        self.plot_win.add_data_point("Uptake", mcs, cell.dict['Uptake'])
        self.plot_win.add_data_point("Secretion", mcs, cell.dict['Secretion'])

        # Do viral model
        for cell in self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING):
            # Step the model for this cell
            CoronavirusLib.step_sbml_model_cell(cell)
            # Pack state variables into cell dictionary
            CoronavirusLib.pack_viral_replication_variables(cell=cell)

            # Test for infection secretion
            if cell.dict['Assembled'] > cell_infection_threshold:
                cell.type = self.INFECTEDSECRETING
                CoronavirusLib.enable_viral_secretion(cell=cell, secretion_rate=secretion_rate)


                # cyttokine params
                cell.dict['ck_production'] = max_ck_secrete_infect  # TODO: replace secretion by hill

            # Test for cell death
            if cell.dict['Assembled'] > cell_death_threshold:
                if not cell.dict['Survived']:
                    p_survival = np.random.random()
                    if p_survival < survival_probability:
                        cell.dict['Survived'] = True
                    else:
                        self.kill_cell(cell=cell)


class Viral_SecretionSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
    """
    def __init__(self, frequency=1):
        CoronavirusSteppableBasePy.__init__(self, frequency)

    def start(self):
        self.track_cell_level_scalar_attribute(field_name='Uptake', attribute_name='Uptake')
        self.track_cell_level_scalar_attribute(field_name='Assembled', attribute_name='Assembled')
        self.track_cell_level_scalar_attribute(field_name='Unpacking', attribute_name='Unpacking')
        self.track_cell_level_scalar_attribute(field_name='Replicating', attribute_name='Replicating')
        self.track_cell_level_scalar_attribute(field_name='Uptake', attribute_name='Uptake')
        self.track_cell_level_scalar_attribute(field_name='Secretion', attribute_name='Secretion')

    def step(self, mcs):
        secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.INFECTEDSECRETING):

            # Evaluate probability of cell uptake of viral particles from environment
            # If cell isn't infected, it changes type to infected here if uptake occurs
            if self.cell_uptakes_virus(viral_field=self.field.Virus,
                                       cell=cell,
                                       diss_coeff_uptake_pr=diss_coeff_uptake_pr,
                                       hill_coeff_uptake_pr=hill_coeff_uptake_pr):
                uptake = secretor.uptakeInsideCellTotalCount(cell,
                                                             cell_infection_threshold / cell.volume,
                                                             relative_viral_uptake)
                cell.dict['Uptake'] = abs(uptake.tot_amount)
                if cell.type == self.UNINFECTED:
                    cell.type = self.INFECTED
                    cell.dict['ck_production'] = max_ck_secrete_infect
                    self.load_viral_replication_model(cell=cell, vr_step_size=vr_step_size,
                                                      unpacking_rate=unpacking_rate,
                                                      replicating_rate=replicating_rate,
                                                      translating_rate=translating_rate,
                                                      packing_rate=packing_rate,
                                                      secretion_rate=secretion_rate)
                CoronavirusLib.set_viral_replication_cell_uptake(cell=cell, uptake=cell.dict['Uptake'])


            if cell.type == self.INFECTEDSECRETING:
                sec_amount = CoronavirusLib.get_viral_replication_cell_secretion(cell=cell)
                secretor.secreteInsideCellTotalCount(cell, sec_amount / cell.volume)


class ImmuneCellKillingSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
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


class ChemotaxisSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
    """
    def __init__(self, frequency=1):
        CoronavirusSteppableBasePy.__init__(self, frequency)

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


class ImmuneCellSeedingSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
    """
    def __init__(self, frequency=1):
        CoronavirusSteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        num_cells = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED, self.INFECTEDSECRETING))
        num_infected = len(self.cell_list_by_type(self.INFECTED, self.INFECTEDSECRETING))
        fraction_infected = 1.0 / num_cells
        if num_cells:
            fraction_infected = num_infected / num_cells

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            p_immune_dying = np.random.random()
            if p_immune_dying < immunecell_dying_rate * fraction_infected:
                cell.targetVolume = 0.0

        if num_infected > 1:
            p_immune_seeding = np.random.random()
            if p_immune_seeding < immune_seeding_rate:
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
                    cell = self.new_cell(self.IMMUNECELL)
                    cell.dict['activated'] = False  # flag for immune cell being naive or activated
                    # cytokine parameters
                    cell.dict['ck_production'] = max_ck_secrete_im  # TODO: replace secretion by hill
                    cell.dict['ck_consumption'] = max_ck_consume  # TODO: replace by hill
                    cell.dict['tot_ck_upt'] = 0

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

        self.pop_data_win = None
        self.pop_data_path = None

        self.med_viral_data_win = None
        self.med_viral_data_path = None

        self.plot_pop_data = plot_pop_data_freq > 0
        self.write_pop_data = write_pop_data_freq > 0

        self.plot_med_viral_data = plot_med_viral_data_freq > 0
        self.write_med_viral_data = write_med_viral_data_freq > 0
        self.med_viral_key = "MedViral"

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
            self.pop_data_win.add_plot("InfectedSecreting", style='Dots', color='green', size=5)
            self.pop_data_win.add_plot("Dying", style='Dots', color='yellow', size=5)
            self.pop_data_win.add_plot("ImmuneCell", style='Dots', color='white', size=5)

        if self.plot_med_viral_data:
            self.med_viral_data_win = self.add_new_plot_window(title='Total diffusive virus',
                                                               x_axis_title='MCS',
                                                               y_axis_title='Number of diffusive viral particles',
                                                               x_scale_type='linear',
                                                               y_scale_type='log',
                                                               grid=True)

            self.med_viral_data_win.add_plot(self.med_viral_key, style='Dots', color='red', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath('pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

            if self.write_med_viral_data:
                self.med_viral_data_path = Path(self.output_dir).joinpath('med_viral_data.dat')
                with open(self.med_viral_data_path, 'w'):
                    pass

    def step(self, mcs):

        plot_pop_data = self.plot_pop_data and mcs % plot_pop_data_freq == 0
        plot_med_viral_data = self.plot_med_viral_data and mcs % plot_med_viral_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % write_pop_data_freq == 0
            write_med_viral_data = self.write_med_viral_data and mcs % write_med_viral_data_freq == 0
        else:
            write_pop_data = False
            write_med_viral_data = False

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.UNINFECTED))
            num_cells_infected = len(self.cell_list_by_type(self.INFECTED))
            num_cells_infectedsecreting = len(self.cell_list_by_type(self.INFECTEDSECRETING))
            num_cells_dying = len(self.cell_list_by_type(self.DYING))
            num_cells_immune = len(self.cell_list_by_type(self.IMMUNECELL))

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

            # Write population data to file if requested
            if write_pop_data:
                with open(self.pop_data_path, 'a') as fout:
                    fout.write('{}, {}, {}, {}, {}, {}\n'.format(mcs,
                                                                 num_cells_uninfected,
                                                                 num_cells_infected,
                                                                 num_cells_infectedsecreting,
                                                                 num_cells_dying,
                                                                 num_cells_immune))

        if plot_med_viral_data or write_med_viral_data:

            # Gather total diffusive viral amount
            med_viral_total = 0.0
            for x, y, z in self.every_pixel():
                med_viral_total += self.field.Virus[x, y, z]

            # Plot total diffusive viral amount if requested
            if plot_med_viral_data and med_viral_total > 0:
                self.med_viral_data_win.add_data_point(self.med_viral_key, mcs, med_viral_total)

            # Write total diffusive viral amount if requested
            if write_med_viral_data:
                with open(self.med_viral_data_path, 'a') as fout:
                    fout.write('{}, {}\n'.format(mcs, med_viral_total))


class CytokineProductionAbsorptionSteppable(CoronavirusSteppableBasePy):
    """
    DESCRIPTION HERE!
    """
    def __init__(self, frequency=1):
        CoronavirusSteppableBasePy.__init__(self, frequency)
        self.track_cell_level_scalar_attribute(field_name='activated', attribute_name='activated')
        self.ck_secretor = None
        self.virus_secretor = None

    def start(self):
        # cytokine diff parameters
        self.get_xml_element('cytokine_dc').cdata = cytokine_dc
        self.get_xml_element('cytokine_decay').cdata = cytokine_field_decay  # no "natural" decay, only consumption
        # and "leakage" outside of simulation laticce

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            # TODO: differentiate this rates between the cell types
            # cytokine production/uptake parameters for immune cells

            cell.dict['ck_production'] = max_ck_secrete_im  # TODO: replace secretion by hill
            cell.dict['ck_consumption'] = max_ck_consume  # TODO: replace by hill

        for cell in self.cell_list_by_type(self.INFECTED,self.INFECTEDSECRETING):
            cell.dict['ck_production'] = max_ck_secrete_infect  # TODO: replace secretion by hill

        # Make sure Secretion plugin is loaded
        # make sure this field is defined in one of the PDE solvers
        # you may reuse secretor for many cells. Simply define it outside the loop
        self.ck_secretor = self.get_field_secretor("cytokine")
        self.virus_secretor = self.get_field_secretor("Virus")

    def step(self, mcs):


        for cell in self.cell_list_by_type(self.INFECTED,self.INFECTEDSECRETING):
            res = self.ck_secretor.secreteInsideCellTotalCount(cell,
                                                               cell.dict['ck_production'] / cell.volume)

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            
            self.virus_secretor.uptakeInsideCellTotalCount(cell,cell.dict['ck_consumption'] / cell.volume, 0.1)
            
            # print(EC50_ck_immune)
            up_res = self.ck_secretor.uptakeInsideCellTotalCount(cell,
                                                                 cell.dict['ck_consumption'] / cell.volume, 0.1)
            #decay seen ck
            cell.dict['tot_ck_upt'] *= ck_memory_immune
            
            #uptake ck
            
            cell.dict['tot_ck_upt'] -= up_res.tot_amount  # from POV of secretion uptake is negative
#             print('tot_upt', cell.dict['tot_ck_upt'],'upt_now', up_res.tot_amount)
            
            p_activate = nCoVUtils.hill_equation(cell.dict['tot_ck_upt'], EC50_ck_immune, 2)
            
            print('prob activation', p_activate, 'upt/ec50', cell.dict['tot_ck_upt']/EC50_ck_immune)
            
            if rng.uniform() < p_activate and not cell.dict['activated'] :
                cell.dict['activated'] = True
                cell.dict['time_activation'] = mcs
            elif (cell.dict['activated'] 
                    and mcs - cell.dict['time_activation'] > minimum_activated_time):
                cell.dict['activated'] = False
                cell.dict['time_activation'] = - 99
            
            
            if cell.dict['activated']:
                # print('activated', cell.id)
                sec_res = self.ck_secretor.secreteInsideCellTotalCount(cell,
                                                                       cell.dict['ck_production'] / cell.volume)
