###############################################################################################################
# To cite this model please use the following:
#
# Josua Aponte-Serrano, T.J. Sego, James A. Glazier,
# "Model of Viral Tissue Infection"
# https://github.com/covid-tissue-models/covid-tissue-response-models/tree/master/tellurium/simple_virus_model
###############################################################################################################

from cc3d.core.PySteppables import *
from cc3d.cpp import CompuCell
import numpy as np

vrl_key = 'viral_replication_loaded'  # Internal use; do not remove

# Conversion Factors
s_to_mcs = 120.0  # s/mcs
um_to_lat_width = 4.0  # um/lattice_length

# Experimental Parameters
exp_cell_diameter = 12.0  # um

exp_replicating_rate = 1.0 / 20.0 * 1.0 / 60.0  # 1.0/20.0min * 1.0min/60.0s = 1.0/1200.0s
exp_translating_rate = exp_replicating_rate * 2.0  #
exp_unpacking_rate = exp_replicating_rate * 20.0
exp_packing_rate = exp_replicating_rate * 4.0
exp_secretion_rate = exp_replicating_rate * 4.0

exp_virus_dc = 10.0 / 100.0  # um^2/s

# cytokines:
# data from https://www.sciencedirect.com/science/article/pii/S1074761317300924 supplemental materials (A)
# and
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3433682/ (B)
cytoplasm_density = 6 # unit = [density of water](B)

exp_cytokine_dc_w = 100 # um^2/s; diffusion constant in water (A,B)
exp_cytokine_dc_cyto = 16 # um^2/s; estimated diffusion constant in cytoplasm (B)

exp_max_cytokine_consumption = 1 # molecule / (cell second); maximum consumption of cytokine; actually a range [0.3,1] molecule / (cell second) (A)
exp_max_cytokine_immune_secretion = 10 # molecule / (cell second) (B)

exp_max_cytokine_consumption_mol = 3.5e-4 # pM/s 
exp_max_cytokine_immune_secretion_mol = 3.5e-3 # pM/s

exp_EC50_cytokine_immune = 50 # pM from (B)


##=============================
# CompuCell Parameters
## cell
cell_diameter = exp_cell_diameter * 1.0 / um_to_lat_width
cell_volume = cell_diameter ** 2
## virus diffusion
virus_dc = exp_virus_dc * s_to_mcs / (um_to_lat_width ** 2)   # virus diffusion constant
virus_dl = cell_diameter * 3.0   # virus diffusion length
virus_decay = virus_dc / (virus_dl ** 2)   # virus decay rate
## virus intra-cellular
unpacking_rate = exp_unpacking_rate * s_to_mcs
replicating_rate = exp_replicating_rate * s_to_mcs
translating_rate = exp_translating_rate * s_to_mcs
packing_rate = exp_packing_rate * s_to_mcs
secretion_rate = exp_secretion_rate * s_to_mcs

## cytokine

cytokine_dc = exp_cytokine_dc_cyto * s_to_mcs / (um_to_lat_width ** 2) # CK diff cst

# pM = pmol/L = pmol/(10^15 um^3) = 10^-15 pmol/(um^3) = 10^-15 * um_to_lat_width^3 pmol/pixel
# pM/s = pM * s_to_mcs / MCS
max_ck_consume = exp_max_cytokine_consumption_mol * um_to_lat_width**3 * s_to_mcs # 1e-15 * pmol/(pixel seconds)
max_ck_secrete_im = exp_max_cytokine_immune_secretion_mol * um_to_lat_width**3 * s_to_mcs # 1e-15 * pmol/(pixel seconds)
EC50_ck_immune = EC50_cytokine_immune * um_to_lat_width**3 # 1e-15 * pmol/pixel

# Threshold at which cell infection is evaluated
cell_infection_threshold = 1.0
# Threshold at which cell death is evaluated
cell_death_threshold = 6.0
# Probability of survival of infected cell once cell_death_threshold is reached
survival_probability = 0.95
# Probability of recovery of infected cell below cell_death_threshold is reached
recovery_probability = 0.001

# Hill equation coefficients for probability of viral particle uptake from the environment
# Measurements are taken w.r.t. the total amount of viral particles in a cell's simulation subdomain
# dissociationt constant
diss_coeff_uptake_pr = 0.5
# Hill coefficient
hill_coeff_uptake_pr = 3.0

# Number of immune cells to seed at the beginning of the simulation
initial_immune_seeding = 10.0
# Rate for seeding of immune cells (constant)
immune_seeding_rate = 1.0 / 10.0
# Max dying rate of immune cells (actual rate is proportional to fraction of infected cells)
immunecell_dying_rate = 1.0 / 500.0

# Name of Antimony/SBML model
vr_model_name = 'viralReplication'
# Antimony/SBML model step size
vr_step_size = 1.0
# Mapping from CellG instance dictionary keys to Antimony/SBML symbols
vr_cell_dict_to_sym = {'Unpacking': 'U',
                       'Replicating': 'R',
                       'Packing': 'P',
                       'Assembled': 'A'}


# Antimony model string generator
# To change models, modify according to this structure
# Modifications should be reflected in
#   1. the items of the dictionary "vr_cell_dict_to_sym"
#   2. the signature of the function "load_viral_replication_model"
# Variable "Uptake" is the uptake variable of the model, and should not be modified
# Variable "Secretion" is the secretion variable of the model, and should not be modified
def viral_replication_model_string(_unpacking_rate, _replicating_rate, _translating_rate, _packing_rate,
                                   _secretion_rate, _u_ini, _r_ini, _p_ini, _a_ini, _uptake=0):
    model_string = """model {}()
      -> U ; Uptake
    U -> R ; unpacking_rate * U;
      -> R ; replicating_rate * R / (0.1 + R);
    R -> P ; translating_rate * R;
    P -> A ; packing_rate * P;
    A -> Secretion ; secretion_rate * A;

    unpacking_rate = {};
    replicating_rate = {};
    translating_rate = {};
    packing_rate = {};
    secretion_rate = {};
    U = {};
    R = {};
    P = {};
    A = {};
    Uptake = {};
    Secretion = 0;
    end""".format(vr_model_name,
                  _unpacking_rate, _replicating_rate, _translating_rate, _packing_rate, _secretion_rate,
                  _u_ini, _r_ini, _p_ini, _a_ini, _uptake)
    return model_string


# Loads viral replication model for a cell
def load_viral_replication_model(_steppable, _cell):
    if _cell.dict[vrl_key]:
        _steppable.delete_sbml_from_cell(vr_model_name, _cell)

    model_string = viral_replication_model_string(unpacking_rate,
                                                  replicating_rate,
                                                  translating_rate,
                                                  packing_rate,
                                                  secretion_rate,
                                                  _cell.dict['Unpacking'],
                                                  _cell.dict['Replicating'],
                                                  _cell.dict['Packing'],
                                                  _cell.dict['Assembled'],
                                                  _cell.dict['Uptake'])
    _steppable.add_antimony_to_cell(model_string=model_string,
                                    model_name=vr_model_name,
                                    cell=_cell,
                                    step_size=vr_step_size)
    _cell.dict[vrl_key] = True
    enable_viral_secretion(_cell, _cell.type == _steppable.INFECTED)


# Calculates the probability of viral uptake from the environment as a function of local viral particle amount
# Returns true if cell uptakes virus
def cell_uptakes_virus(_steppable, viral_field, _cell):
    # Calculate total viral amount in cell's domain
    cell_env_viral_val = 0.0
    for ptd in _steppable.get_cell_pixel_list(_cell):
        cell_env_viral_val += viral_field[ptd.pixel.x, ptd.pixel.y, ptd.pixel.z]

    # Evaluate probability of uptake
    if cell_env_viral_val != 0:
        max_uptake_pr = 1 / (1 + (diss_coeff_uptake_pr / cell_env_viral_val) ** hill_coeff_uptake_pr)
        return np.random.random() < max_uptake_pr
    else:
        return False


# These are prototypes of specialized utility functions for this project; do not modify

# Removes viral replication model for a cell
def remove_viral_replication_model(_steppable, _cell):
    if _cell.dict[vrl_key]:
        _steppable.delete_sbml_from_cell(vr_model_name, _cell)
        _cell.dict[vrl_key] = False


# Sets the current state variable "Uptake" for a cell
def set_viral_replication_cell_uptake(_cell, _uptake):
    assert _cell.dict[vrl_key]
    getattr(_cell.sbml, vr_model_name)['Uptake'] = _uptake


# Gets the current state variable "Secretion" for a cell
def get_viral_replication_cell_secretion(_cell):
    assert _cell.dict[vrl_key]
    secr = getattr(_cell.sbml, vr_model_name)['Secretion']
    getattr(_cell.sbml, vr_model_name)['Secretion'] = 0.0
    return secr


# Enable/disable secretion in intracellular model
def enable_viral_secretion(_cell, _enable: bool = True):
    if _enable:
        getattr(_cell.sbml, vr_model_name)['secretion_rate'] = secretion_rate
    else:
        getattr(_cell.sbml, vr_model_name)['secretion_rate'] = 0.0


# Loads state variables from SBML into cell dictionary
def pack_viral_replication_variables(_cell):
    assert _cell.dict[vrl_key]
    for k, v in vr_cell_dict_to_sym.items():
        _cell.dict[k] = getattr(_cell.sbml, vr_model_name)[v]


# Loads state variables from SBML into cell dictionary
def reset_viral_replication_variables(_cell):
    _cell.dict['Uptake'] = 0
    _cell.dict['Secretion'] = 0
    for k in vr_cell_dict_to_sym.keys():
        _cell.dict[k] = 0


# Steps SBML model for a cell
def step_viral_replication_cell(_cell):
    dict_attrib = CompuCell.getPyAttrib(_cell)
    assert 'SBMLSolver' in dict_attrib
    dict_attrib['SBMLSolver'][vr_model_name].timestep()


class CellsInitializerSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.get_xml_element('virus_dc').cdata = virus_dc
        self.get_xml_element('virus_decay').cdata = virus_decay
        
        # cytokine diff parameters        
        self.get_xml_element('cytokine_dc').cdata = cytokine_dc
        self.get_xml_element('cytokine_decay').cdata = 0 # no "natural" decay, only consumption
        
        
        for x in range(0, self.dim.x, int(cell_diameter)):
            for y in range(0, self.dim.y, int(cell_diameter)):
                cell = self.new_cell(self.UNINFECTED)
                self.cellField[x:x + int(cell_diameter), y:y + int(cell_diameter), 0] = cell
                cell.dict[vrl_key] = False
                reset_viral_replication_variables(cell)
                cell.dict['Survived'] = False
                cell.dict['Recovered'] = False
                if x == int(self.dim.x / 2):
                    if y == int(self.dim.x / 2):
                        # Start infection of an uninfected cell:
                        cell.dict['Unpacking'] = 1.0

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
            # cytokine production/uptake parameters for immune cells
            cell.dict['immune_production'] = max_ck_secrete_im ##TODO: replace secretion by hill
            cell.dict['immune_consumption'] = max_ck_consume ##TODO: replace by hill


class Viral_ReplicationSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

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
        for cell in self.cell_list_by_type(self.UNINFECTED, self.INFECTED):
            load_viral_replication_model(self, cell)

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
        for cell in self.cell_list_by_type(self.UNINFECTED, self.INFECTED):
            # Step the model for this cell
            step_viral_replication_cell(cell)
            # Pack state variables into cell dictionary
            pack_viral_replication_variables(cell)

            # Test for infection
            if cell.dict['Assembled'] > cell_infection_threshold:
                cell.type = self.INFECTED
                enable_viral_secretion(cell)

            # Test for cell death
            if cell.dict['Assembled'] > cell_death_threshold:
                if not cell.dict['Survived']:
                    p_survival = np.random.random()
                    if p_survival < survival_probability:
                        cell.dict['Survived'] = True
                    else:
                        cell.type = self.DYING
                        reset_viral_replication_variables(cell)
                        # Remove state model: no model for dead cell type
                        remove_viral_replication_model(self, cell)


class Viral_SecretionSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        self.track_cell_level_scalar_attribute(field_name='Uptake', attribute_name='Uptake')
        self.track_cell_level_scalar_attribute(field_name='Assembled', attribute_name='Assembled')
        self.track_cell_level_scalar_attribute(field_name='Unpacking', attribute_name='Unpacking')
        self.track_cell_level_scalar_attribute(field_name='Replicating', attribute_name='Replicating')
        self.track_cell_level_scalar_attribute(field_name='Uptake', attribute_name='Uptake')
        self.track_cell_level_scalar_attribute(field_name='Secretion', attribute_name='Secretion')

    def step(self, mcs):
        secretor = self.get_field_secretor("Virus")
        for cell in self.cell_list_by_type(self.INFECTED, self.UNINFECTED):
            relative_viral_uptake = 0.1

            # Evaluate probability of cell uptake of viral particles from environment
            if cell_uptakes_virus(self, self.field.Virus, cell):
                uptake = secretor.uptakeInsideCellTotalCount(cell,
                                                             cell_infection_threshold / cell.volume,
                                                             relative_viral_uptake)
                cell.dict['Uptake'] = abs(uptake.tot_amount)
                set_viral_replication_cell_uptake(cell, cell.dict['Uptake'])

        for cell in self.cell_list_by_type(self.INFECTED):
            sec_amount = get_viral_replication_cell_secretion(cell)
            secretor.secreteInsideCellTotalCount(cell, sec_amount / cell.volume)


class ImmuneCellKillingSteppable(SteppableBasePy):
    def step(self, mcs):
        for cell in self.cell_list_by_type(self.INFECTED):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.IMMUNECELL:
                        cell.type = self.DYING
                        reset_viral_replication_variables(cell)
                        # Remove state model
                        remove_viral_replication_model(self, cell)


class ChemotaxisSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        for cell in self.cell_list_by_type(self.IMMUNECELL):
            cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
            cd.setLambda(50.0)
            cd.assignChemotactTowardsVectorTypes([self.MEDIUM])

    def step(self, mcs):
        field = self.field.Virus
        for cell in self.cell_list_by_type(self.IMMUNECELL):
            cd = self.chemotaxisPlugin.getChemotaxisData(cell, "Virus")
            concentration = field[cell.xCOM, cell.yCOM, 0]
            constant = 50.0
            l = constant / (1.0 + concentration)
            cd.setLambda(l)


class RecoverySteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        for cell in self.cell_list_by_type(self.INFECTED):
            if not cell.dict['Recovered'] and not cell.dict['Survived']:
                p_recovery = np.random.random()
                if p_recovery < recovery_probability:
                    cell.type = self.UNINFECTED
                    reset_viral_replication_variables(cell)
                    cell.dict['Recovered'] = True
                    # Update state model
                    load_viral_replication_model(self, cell)


class ImmuneCellSeedingSteppable(SteppableBasePy):
    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def step(self, mcs):
        num_cells = len(self.cell_list_by_type(self.UNINFECTED, self.INFECTED))
        num_infected = len(self.cell_list_by_type(self.INFECTED))
        fraction_infected = 1.0 / num_cells
        if num_cells:
            fraction_infected = num_infected / num_cells

        for cell in self.cell_list_by_type(self.IMMUNECELL):
            p_immune_dying = np.random.random()
            if p_immune_dying < immunecell_dying_rate * fraction_infected:
                cell.targetVolume = 0.0

        p_immune_seeding = np.random.random()
        if p_immune_seeding < immune_seeding_rate:
            open_space = True
            viral_concentration = 0
            for iteration in range(10):
                xi = np.random.randint(0, self.dim.x - 2 * cell_diameter)
                yi = np.random.randint(0, self.dim.y - 2 * cell_diameter)
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
                self.cellField[x_seed:x_seed + int(cell_diameter), y_seed:y_seed + int(cell_diameter), 1] = cell
                cd = self.chemotaxisPlugin.addChemotaxisData(cell, "Virus")
                cd.setLambda(50.0)
                cd.assignChemotactTowardsVectorTypes([self.MEDIUM])
                cell.targetVolume = cell_volume
                cell.lambdaVolume = cell_volume
