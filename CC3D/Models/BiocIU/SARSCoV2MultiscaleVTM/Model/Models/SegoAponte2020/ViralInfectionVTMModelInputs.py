# Info for documentation of parameter values goes here, in *__param_desc__*
# The parameter name and assigned description can be quickly exported to csv along with specified values using
# *export_parameters*, in export_parameters.py
# To use, assign the name of the code variable as a key, and a string description as its value
# e.g.,     __param_desc__['my_var'] = 'My variable'
#           my_var = 1
# The generated csv will read, "my_var, 1, My variable'
__param_desc__ = {}

# Data control options
__param_desc__['track_model_variables'] = 'Enables cell-level tracking of model variables (for rendering)'
track_model_variables = False  # Set to true to enable cell-level tracking of model variables (for rendering)
__param_desc__['plot_pop_data_freq'] = 'Plot population data frequency'
plot_pop_data_freq = 0  # Plot population data frequency (disable with 0)
__param_desc__['write_pop_data_freq'] = 'Write population data to simulation directory frequency'
write_pop_data_freq = 0  # Write population data to simulation directory frequency (disable with 0)
__param_desc__['plot_med_diff_data_freq'] = 'Plot total diffusive field amount frequency'
plot_med_diff_data_freq = 0  # Plot total diffusive field amount frequency (disable with 0)
__param_desc__['write_med_diff_data_freq'] = 'Write total diffusive field amount frequency'
write_med_diff_data_freq = 0  # Write total diffusive field amount frequency (disable with 0)
__param_desc__['plot_ir_data_freq'] = 'Plot immune recruitment data frequency'
plot_ir_data_freq = 0  # Plot immune recruitment data frequency (disable with 0)
__param_desc__['write_ir_data_freq'] = 'Write immune recruitment data to simulation directory frequency'
write_ir_data_freq = 0  # Write immune recruitment data to simulation directory frequency (disable with 0)
__param_desc__['plot_death_data_freq'] = 'Plot death data frequency'
plot_death_data_freq = 0  # Plot death data frequency (disable with 0)
__param_desc__['write_death_data_freq'] = 'Write death data to simulation directory frequency'
write_death_data_freq = 0  # Write death data to simulation directory frequency (disable with 0)

# Helpful conversions

# pM = pmol/L = pmol/(10^15 um^3) = 10^-15 pmol/(um^3) = 10^-15 * voxel_length^3 pmol/pixel
# pM/s = pM * step_period / MCS

__param_desc__['mol_p_L_2_mol_p_um3'] = 'Conversion from mol/L to mol/um^3'
mol_p_L_2_mol_p_um3 = 1e-15

__param_desc__['pmol_to_cc3d_au'] = 'Scale factor for concentration'
pmol_to_cc3d_au = 1e14  # 1e14au/1pmol

# Experimental Parameters
__param_desc__['cell_diameter'] = 'Cell diameter'
cell_diameter = 12.0  # um
__param_desc__['cell_volume'] = 'Cell volume'
cell_volume = cell_diameter ** 2  # um2

__param_desc__['unpacking_rate'] = 'Unpacking rate'
unpacking_rate = 1.0 / 100.0 * 1.0 / 60.0
__param_desc__['replicating_rate'] = 'Replicating rate'
replicating_rate = 1.0 / 200.0 * 1.0 / 60.0  # 1.0/20.0min * 1.0min/60.0s = 1.0/1200.0s
__param_desc__['translating_rate'] = 'Translating rate'
translating_rate = 1.0 / 300.0 * 1.0 / 60.0
__param_desc__['packing_rate'] = 'Packing rate'
packing_rate = 1.0 / 100.0 * 1.0 / 60.0
__param_desc__['secretion_rate'] = 'Secretion rate of assembled virus'
secretion_rate = 1.0 / 100.0 * 1.0 / 60.0

__param_desc__['virus_dc'] = 'Viral diffusion coefficient'
virus_dc = 10.0 / 1000.0  # um^2/s
__param_desc__['virus_dl'] = 'Unitless virus diffusion length'
virus_dl = cell_diameter * 3.0  # virus diffusion length
__param_desc__['virus_decay'] = 'Unitless virus decay coefficient'
virus_decay = virus_dc / (virus_dl ** 2)  # virus decay rate

# cytokines:
# data from https://www.sciencedirect.com/science/article/pii/S1074761317300924 supplemental materials (A) and
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3433682/ (B)
# cytoplasm_density = 6  # unit = [density of water](B)
# 100  um^2/s; diffusion constant in water (A,B)

__param_desc__['cytokine_dc'] = 'Cytokine diffusion coefficient'
cytokine_dc = 16 / 100  # um^2/s; estimated diffusion constant in cytoplasm (B)

__param_desc__['max_ck_diff_len'] = 'Cytokine diffusion length'
max_ck_diff_len = 100  # um  from (A); this diffusion length is due to uptake, I'm using it as a base
# for the decay due to ck leakage outside of the simulated lattice.
# dl = sqrt(D/g) -> g = D/dl**2
__param_desc__['cytokine_field_decay'] = 'Cytokine decay coefficient'
cytokine_field_decay = cytokine_dc/(1.1*max_ck_diff_len)**2

# molecule / (cell second); maximum consumption of cytokine; actually a range [0.3,1] molecule / (cell second) (A)
__param_desc__['max_ck_consume'] = 'Cytokine immune uptake rate'
max_ck_consume = 3.5e-4  # pM/s (B) -- they also have it in 1 molecule / (cell second)
__param_desc__['max_ck_secrete_im'] = 'Maximum cytokine immune cell secretion rate'
max_ck_secrete_im = 3.5e-4  # pM/s (B) -- they also have it in 10 molecule / (cell second)

__param_desc__['EC50_ck_immune'] = 'Immune cell cytokine activation'
EC50_ck_immune = 10  # pM from (B), it's a range from [1,50]pM

__param_desc__['minimum_activated_time'] = 'Immune cell activated time'
minimum_activated_time = 600 * 60  # min * s/min

__param_desc__['exp_max_amount_viral_mRNA'] = 'Scale factor for number of mRNA per infected cell'
exp_max_amount_viral_mRNA = 1000

# 10^5 ~ 10^7 /mL RNA
# ~ 2000 um^3 cell volume -> 2e-9 mL/cell
#

# oxidation agent

__param_desc__['oxi_dl'] = 'Oxidation Agent diffusion length'
oxi_dl = 3 * cell_diameter  # um

__param_desc__['oxi_dc'] = 'Oxidation Agent diffusion coefficient'
oxi_dc = 4 * cytokine_dc  # um^2/s

__param_desc__['oxi_decay'] = 'Oxidation Agent decay coefficient'
oxi_decay = oxi_dc / oxi_dl**2  # 1/s


# =============================
# CompuCell Parameters
# cell
__param_desc__['volume_lm'] = 'Volume constraint LM'
volume_lm = 9

# cytokine / immune activation

__param_desc__['ec50_immune_ck_prod'] = 'Amount of seen cytokine at 50% cytokine production by immune cells'
ec50_immune_ck_prod = 1.0  # pM

__param_desc__['ck_equilibrium'] = 'equilibrium amount of ck in immune surface'
ck_equilibrium = 2.1*EC50_ck_immune  # equilibrium amount of ck in immune surface
__param_desc__['ck_memory_immune'] = '1 - Immune cell bound cytokine memory'
ck_memory_immune = 1 - max_ck_consume/ck_equilibrium  # decay therm for "seen" ck by immune

__param_desc__['max_ck_secrete_infect'] = 'Maximum cytokine lung tissue secretion rate'
max_ck_secrete_infect = 10 * max_ck_secrete_im  # pM/s

__param_desc__['ec50_infecte_ck_prod'] = 'Amount of internal assembled virus to be at 50% cytokine production'
ec50_infecte_ck_prod = 0.1  # amount of 'internal assembled virus' to be at 50% ck production; chosen from
# typical simulation values of cell.dict['Uptake'] + cell.dict['Assembled']. they stay around .1 and go up as the
# simulation progresses

# oxidation agent

__param_desc__['oxi_secr_thr'] = 'Immune cell cytokine concentration threshold for Oxidation Agent release'
oxi_secr_thr = 1.5625  # pM; -> 10 a.u. at 4 um/lattice length

__param_desc__['max_oxi_secrete'] = 'Immune cell oxidation agent secretion rate'
max_oxi_secrete = max_ck_secrete_infect  # pM/s

__param_desc__['oxi_death_thr'] = 'Tissue cell Oxidation Agent threshold for death'
oxi_death_thr = 0.234375  # pM; -> 1.5 a.u. at 4 um/lattice length

# Threshold at which cell infection is evaluated
__param_desc__['cell_infection_threshold'] = 'Threshold of assembled viral particles above which infected become ' \
                                             'infectedSecreting '
cell_infection_threshold = 1.0

# Hill equations coefficients for probability of viral-induced apoptosis
# Measurements are taken w.r.t. the total amount of assembled viral particles in a cell's simulation subdomain
# dissociationt constant
__param_desc__['diss_coeff_uptake_apo'] = 'Dissociation coefficient for probability of viral-induced apoptosis'
diss_coeff_uptake_apo = 100.0
# Hill coefficient
__param_desc__['hill_coeff_uptake_apo'] = 'Hill coefficient for probability of viral-induced apoptosis'
hill_coeff_uptake_apo = 2.0

# Hill equation coefficients for probability of viral particle uptake from the environment
# Measurements are taken w.r.t. the total amount of viral particles in a cell's simulation subdomain
# Hill coefficient
__param_desc__['hill_coeff_uptake_pr'] = 'Hill coefficient for probability of viral particle uptake from the ' \
                                         'environment '
hill_coeff_uptake_pr = 2.0
# Rate coefficient
__param_desc__['rate_coeff_uptake_pr'] = 'Rate coefficient for probability of viral particle uptake from the ' \
                                         'environment '
rate_coeff_uptake_pr = 1200.0  # s

# Efficiency of viral uptake
__param_desc__['relative_viral_uptake'] = 'Efficiency of viral uptake'
relative_viral_uptake = 0.1

# Number of immune cells to seed at the beginning of the simulation
# Leave as zero if running through CallableCC3D (bug in NeighborTracker)
__param_desc__['initial_immune_seeding'] = 'Number of immune cells to seed at the beginning of the simulation'
initial_immune_seeding = 0.0
# Bystander effect
__param_desc__['bystander_effect'] = 'Probability of death due to the Bystander Effect'
bystander_effect = 0.5
# Lambda Chemotaxis
__param_desc__['lamda_chemotaxis'] = 'Lambda chemotaxis (chemotactic sensitivity)'
lamda_chemotaxis = 1.0

# Antimony/SBML model step size
__param_desc__['vr_step_size'] = 'Antimony/SBML model step size'
vr_step_size = 1.0

# Viral Internalization parameters
__param_desc__['kon'] = 'Virus-receptors association affinity'
kon = 1.4E9  # 1/(M * s)
__param_desc__['koff'] = 'Virus-receptors disassociation affinity'
koff = 1.4E-4  # 1/s

# Initial number of receptors on the cell surface
__param_desc__['initial_unbound_receptors'] = 'Initial number of receptors on the cell surface'
initial_unbound_receptors = 200

# Number of cell receptors at which the replication rate is half max (sets the steady state value of mRNA)
__param_desc__['r_half'] = 'Number of cell receptors at which the replication rate is half of maximum'
r_half = exp_max_amount_viral_mRNA/(replicating_rate/translating_rate-1)

# Imunne recruitment parameters

# State variable rate addition coefficient
__param_desc__['ir_add_coeff'] = 'Immune response rate addition coefficient'
ir_add_coeff = 1.0 / 1200.0  # 1/s
# State variable rate subtraction coefficient
__param_desc__['ir_subtract_coeff'] = 'Immune response rate subtraction coefficient'
if initial_immune_seeding == 0:
    ir_subtract_coeff = ir_add_coeff / 5.0 / 1200.0  # 1/cell/s
else:
    ir_subtract_coeff = ir_add_coeff / initial_immune_seeding / 1200.0  # 1/cell/s
# State variable rate delay coefficient
__param_desc__['ir_delay_coeff'] = 'Immune response rate delay coefficient'
ir_delay_coeff = 1E3 * 1200.0  # s * cytokine
# State variable rate decay coefficient
__param_desc__['ir_decay_coeff'] = 'Immune response rate decay coefficient'
ir_decay_coeff = 1E-1 / 1200.0  # 1/s
# Ratio of decay cytokine used as recruitment signal
__param_desc__['ir_transmission_coeff'] = 'Immune response cytokine transmission coefficient'
ir_transmission_coeff = 5E-1  # ul
# Scales state variable in probability functions
__param_desc__['ir_prob_scaling_factor'] = 'Immune response probability scaling coefficient'
ir_prob_scaling_factor = 1.0 / 100.0  # ul
