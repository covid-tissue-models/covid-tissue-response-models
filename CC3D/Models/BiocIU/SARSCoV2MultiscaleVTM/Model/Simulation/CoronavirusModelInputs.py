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
__param_desc__['plot_vrm_data_freq'] = 'Plot viral replication model data frequency'
plot_vrm_data_freq = 0  # Plot viral replication model data frequency (disable with 0)
__param_desc__['write_vrm_data_freq'] = 'Write viral replication model data to simulation directory frequency'
write_vrm_data_freq = 0  # Write viral replication model data to simulation directory frequency (disable with 0)
__param_desc__['plot_vim_data_freq'] = 'Plot viral internalization model data frequency'
plot_vim_data_freq = 0  # Plot viral internalization model data frequency (disable with 0)
__param_desc__['write_vim_data_freq'] = 'Write viral internalization model data to simulation directory frequency'
write_vim_data_freq = 0  # Write viral internalization model data to simulation directory frequency (disable with 0)
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
__param_desc__['plot_spat_data_freq'] = 'Plot spatial data frequency'
plot_spat_data_freq = 0  # Plot spatial data frequency (disable with 0)
__param_desc__['write_spat_data_freq'] = 'Write spatial data to simulation directory frequency'
write_spat_data_freq = 0  # Write spatial data to simulation directory frequency (disable with 0)

# Conversion Factors
__param_desc__['s_to_mcs'] = 'Simulation step'
s_to_mcs = 1200.0  # s/mcs
__param_desc__['um_to_lat_width'] = 'Lattice width'
um_to_lat_width = 4.0  # um/lattice_length

__param_desc__['pmol_to_cc3d_au'] = 'Scale factor for concentration'
pmol_to_cc3d_au = 1e14  # 1e15au/1pmol

# Experimental Parameters
__param_desc__['exp_cell_diameter'] = 'Cell diameter'
exp_cell_diameter = 12.0  # um

__param_desc__['exp_unpacking_rate'] = 'Unpacking rate'
exp_unpacking_rate = 1.0 / 100.0 * 1.0 / 60.0
__param_desc__['exp_replicating_rate'] = 'Replicating rate'
exp_replicating_rate = 1.0 / 200.0 * 1.0 / 60.0  # 1.0/20.0min * 1.0min/60.0s = 1.0/1200.0s
__param_desc__['exp_translating_rate'] = 'Translating rate'
exp_translating_rate = 1.0 / 300.0 * 1.0 / 60.0
__param_desc__['exp_packing_rate'] = 'Packing rate'
exp_packing_rate = 1.0 / 100.0 * 1.0 / 60.0
__param_desc__['exp_secretion_rate'] = 'Secretion rate'
exp_secretion_rate = 1.0 / 100.0 * 1.0 / 60.0

__param_desc__['exp_virus_dc'] = 'Viral diffusion coefficient'
exp_virus_dc = 10.0 / 1000.0  # um^2/s

# cytokines:
# data from https://www.sciencedirect.com/science/article/pii/S1074761317300924 supplemental materials (A)
# and
# from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3433682/ (B)
# cytoplasm_density = 6  # unit = [density of water](B)

# __param_desc__['exp_cytokine_dc_w'] = 'Cytokine extracellular diffusion coefficient'
# exp_cytokine_dc_w = 100  # um^2/s; diffusion constant in water (A,B)
# the division by 10 is due to the small lattice size, as fiscussed with josh. we all should discuss this matter
__param_desc__['exp_cytokine_dc_cyto'] = 'Cytokine diffusion coefficient'
exp_cytokine_dc_cyto = 16 / 100  # um^2/s; estimated diffusion constant in cytoplasm (B)
# ^ the  /10 is not experimental; added because of relative small area and because virus D is (or was) slowed down

__param_desc__['exp_max_ck_diff_len'] = 'Cytokine diffusion length'
exp_max_ck_diff_len = 100  # um  from (A); this diffusion length is due to uptake, I'm using it as a base
# for the decay due to ck leakage outside of the simulatted lattice.
# dl = sqrt(D/g) -> g = D/dl**2
__param_desc__['exp_min_ck_decay'] = 'Cytokine decay coefficient'
exp_min_ck_decay = exp_cytokine_dc_cyto/(1.1*exp_max_ck_diff_len)**2

# molecule / (cell second); maximum consumption of cytokine; actually a range [0.3,1] molecule / (cell second) (A)
__param_desc__['exp_max_cytokine_consumption'] = 'Maximum cytokine consumption'
exp_max_cytokine_consumption = 1
# molecule / (cell second) (B)
__param_desc__['exp_max_cytokine_immune_secretion'] = 'Maximum cytokine immune secretion rate'
exp_max_cytokine_immune_secretion = 10

__param_desc__['exp_max_cytokine_consumption_mol'] = 'Cytokine immune uptake rate'
exp_max_cytokine_consumption_mol = 3.5e-4  # pM/s
__param_desc__['exp_max_cytokine_immune_secretion_mol'] = 'Maximum cytokine lung tissue secretion rate'
exp_max_cytokine_immune_secretion_mol = 3.5e-3  # pM/s

__param_desc__['exp_EC50_cytokine_immune'] = 'Immune cell cytokine activation'
exp_EC50_cytokine_immune = 10  # pM from (B), it's a range from [1,50]pM
# tbd: try to find experimental data
__param_desc__['minimum_activated_time_seconds'] = 'Immune cell activated time'
minimum_activated_time_seconds = 600 * 60  # min * s/min

__param_desc__['exp_max_amount_viral_mRNA'] = 'Scale factor for number of mRNA per infected cell'
exp_max_amount_viral_mRNA = 1000

# oxidation agent

__param_desc__['exp_oxi_dl'] = 'Oxidation Agent diffusion length'
exp_oxi_dl = 3 * exp_cell_diameter  # [um]; guestimation; [.3,3]
# __param_desc__['exp_oxi_dc_water'] = 'Oxidation Agent extracellular diffusion coefficient'
# exp_oxi_dc_water = 1.2  # cm2/day; http://www.idc-online.com/technical_references/pdfs/chemical_engineering/Transport_Properties_of_Hydrogen_Peroxide_II.pdf
# exp_oxi_dc_water = exp_oxi_dc_water * 1e8 / 86400  # um2/s
# exp_oxi_dc_cyto = exp_oxi_dc_water * .16  # rescale by relative density; cyto ~ 6*water
# exp_oxi_dc_cyto = exp_cytokine_dc_cyto

# the experimental values are WAY WAY too high for the simulation to behave properlly, so:
__param_desc__['exp_oxi_dc_cyto'] = 'Oxidation Agent diffusion coefficient'
exp_oxi_dc_cyto = 4 * exp_cytokine_dc_cyto


# =============================
# CompuCell Parameters
# cell
__param_desc__['cell_diameter'] = 'Unitless cell diameter'
cell_diameter = exp_cell_diameter * 1.0 / um_to_lat_width
__param_desc__['cell_volume'] = 'Unitless cell volume'
cell_volume = cell_diameter ** 2
# virus diffusion
__param_desc__['virus_dc'] = 'Unitless virus diffusion coefficient'
virus_dc = exp_virus_dc * s_to_mcs / (um_to_lat_width ** 2)  # virus diffusion constant
__param_desc__['virus_dl'] = 'Unitless virus diffusion length'
virus_dl = cell_diameter * 3.0  # virus diffusion length
__param_desc__['virus_decay'] = 'Unitless virus decay coefficient'
virus_decay = virus_dc / (virus_dl ** 2)  # virus decay rate
# virus intra-cellular
__param_desc__['unpacking_rate'] = 'Unitless unpacking rate'
unpacking_rate = exp_unpacking_rate * s_to_mcs
__param_desc__['replicating_rate'] = 'Unitless replicating rate'
replicating_rate = exp_replicating_rate * s_to_mcs
__param_desc__['translating_rate'] = 'Unitless translating rate'
translating_rate = exp_translating_rate * s_to_mcs
__param_desc__['packing_rate'] = 'Unitless packing rate'
packing_rate = exp_packing_rate * s_to_mcs
__param_desc__['secretion_rate'] = 'Unitless secretion rate'
secretion_rate = exp_secretion_rate * s_to_mcs

# cytokine / immune activation

__param_desc__['cytokine_dc'] = 'Unitless cytokine diffusion coefficient'
cytokine_dc = exp_cytokine_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)  # CK diff cst
# [1/g] = [s] -> [1/g]/s_to_mcs = [s]/[s/mcs] = [mcs]
# -> [g] = [1/s] -> [g]*[s_to_mcs] = [1/s]*[s/mcs] = [1/mcs]
__param_desc__['cytokine_field_decay'] = 'Unitless cytokine decay coefficient'
cytokine_field_decay = exp_min_ck_decay * s_to_mcs
# pM = pmol/L = pmol/(10^15 um^3) = 10^-15 pmol/(um^3) = 10^-15 * um_to_lat_width^3 pmol/pixel
# pM/s = pM * s_to_mcs / MCS
__param_desc__['max_ck_consume'] = 'Unitless maximum cytokine consumption'
max_ck_consume = exp_max_cytokine_consumption_mol * um_to_lat_width ** 3 * s_to_mcs * 1e-15 * pmol_to_cc3d_au  # cc3d_au/(pixel seconds)
__param_desc__['max_ck_secrete_im'] = 'Unitless maximum cytokine immune secretion rate'
max_ck_secrete_im = exp_max_cytokine_immune_secretion_mol * um_to_lat_width ** 3 * s_to_mcs * 1e-15 * pmol_to_cc3d_au  # * cc3d_au/(pixel seconds)
__param_desc__['EC50_ck_immune'] = 'Unitless immune cell cytokine activation'
EC50_ck_immune = exp_EC50_cytokine_immune * um_to_lat_width ** 3 * 1e-15 * pmol_to_cc3d_au  # * cc3d_au/pixel
__param_desc__['ck_equilibrium'] = 'equilibrium amount of ck in immune surface'
# ck_equilibrium = 1.5*EC50_ck_immune # equilibrium amount of ck in immune surface
ck_equilibrium = 2.1*EC50_ck_immune  # equilibrium amount of ck in immune surface
__param_desc__['ck_memory_immune'] = '1 - Immune cell bound cytokine memory'
ck_memory_immune = 1 - max_ck_consume/ck_equilibrium  # decay therm for "seen" ck by immune

__param_desc__['max_ck_secrete_infect'] = 'Unitless maximum cytokine lung tissue secretion rate'
max_ck_secrete_infect = 10*max_ck_secrete_im
__param_desc__['minimum_activated_time'] = 'Unitless immune cell activated time'
minimum_activated_time = minimum_activated_time_seconds/s_to_mcs  # mcs

__param_desc__['ec50_infecte_ck_prod'] = 'Amount of internal assembled virus to be at 50% cytokine production'
ec50_infecte_ck_prod = 0.1  # amount of 'internal assembled virus' to be at 50% ck production; chosen from
# tipical simulation values of cell.dict['Uptake'] + cell.dict['Assembled']. they stay around .1 and go up as the
# simulation progresses

# oxidation agent


__param_desc__['oxi_dl'] = 'Unitless Oxidation Agent diffusion length'
oxi_dl = exp_oxi_dl/um_to_lat_width

__param_desc__['oxi_dc'] = 'Unitless Oxidation Agent diffusion coefficient'
oxi_dc = exp_oxi_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)

__param_desc__['oxi_decay'] = 'Unitless Oxidation Agent decay coefficient'
oxi_decay = oxi_dl**2/oxi_dc

__param_desc__['oxi_sec_thr'] = 'Immune cell cytokine concentration threshold for Oxidation Agent release'
oxi_sec_thr = 10

__param_desc__['max_oxi_secrete'] = 'Immune cell oxidation agent secretion rate'
max_oxi_secrete = max_ck_secrete_infect

__param_desc__['oxi_death_thr'] = 'Tissue cell Oxidation Agent threshold for death'
oxi_death_thr = 1.5


# Threshold at which cell infection is evaluated
__param_desc__['cell_infection_threshold'] = 'Threshold of assembled viral particles above which infected become infectedSecreting'
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
__param_desc__['hill_coeff_uptake_pr'] = 'Hill coefficient for probability of viral particle uptake from the environment'
hill_coeff_uptake_pr = 2.0
# Rate coefficient
__param_desc__['rate_coeff_uptake_pr'] = 'Rate coefficient for probability of viral particle uptake from the environment'
rate_coeff_uptake_pr = 1200.0

# Efficiency of viral uptake
__param_desc__['relative_viral_uptake'] = 'Efficiency of viral uptake'
relative_viral_uptake = 0.1

# Number of immune cells to seed at the beginning of the simulation
# Leave as zero if running through CallableCC3D (bug in NeighborTracker)
__param_desc__['initial_immune_seeding'] = 'Number of immune cells to seed at the beginning of the simulation'
initial_immune_seeding = 0.0
# Bystander effect
__param_desc__['bystander_effect'] = 'Probability rate of death due to the Bystander Effect per MCS'
bystander_effect = 2.0/4.0
# Lambda Chemotaxis
__param_desc__['lamda_chemotaxis'] = 'Lambda chemotaxis (chemotactic sensitivity)'
# lamda_chemotaxis = 100.0
lamda_chemotaxis = 100.0/100.0

# Antimony/SBML model step size
__param_desc__['vr_step_size'] = 'Antimony/SBML model step size'
vr_step_size = 1.0

# Viral Internalization parameters
__param_desc__['exp_kon'] = 'Virus-receptors association affinity'
exp_kon = 1.4E4  # 1/(M * s)
__param_desc__['exp_koff'] = 'Virus-receptors disassociation affinity'
exp_koff = 1.4E-4  # 1/s

# Initial number of receptors on the cell surface
__param_desc__['initial_unbound_receptors'] = 'Initial number of receptors on the cell surface'
initial_unbound_receptors = 200

# Number of cell receptors at which the replication rate is half max (sets the steady state value of mRNA)
__param_desc__['r_half'] = 'Number of cell receptors at which the replication rate is half of maximum'
r_half = exp_max_amount_viral_mRNA/(replicating_rate/translating_rate-1)

__param_desc__['kon'] = 'Unitless virus-receptors association affinity'
kon = exp_kon * s_to_mcs * 1.0E15 * (1.0/(um_to_lat_width**3)) * (1.0/1.0E12) * (1.0/pmol_to_cc3d_au) * 100
__param_desc__['koff'] = 'Unitless virus-receptors disassociation affinity'
koff = exp_koff * s_to_mcs

# Imunne recruitment parameters

# State variable rate addition coefficient
__param_desc__['ir_add_coeff'] = 'Immune response rate addition coefficient'
ir_add_coeff = 1.0
# State variable rate subtraction coefficient
__param_desc__['ir_subtract_coeff'] = 'Immune response rate subtraction coefficient'
if initial_immune_seeding == 0:
    ir_subtract_coeff = ir_add_coeff / 5.0
else:
    ir_subtract_coeff = ir_add_coeff / initial_immune_seeding
# State variable rate delay coefficient
__param_desc__['ir_delay_coeff'] = 'Immune response rate delay coefficient'
ir_delay_coeff = 1*1E3
# State variable rate decay coefficient
__param_desc__['ir_decay_coeff'] = 'Immune response rate decay coefficient'
ir_decay_coeff = 1E-1
# Ratio of decay cytokine used as recruitment signal
__param_desc__['ir_transmission_coeff'] = 'Immune response cytokine transmission coefficient'
ir_transmission_coeff = 5E-1
# Scales state variable in probability functions
__param_desc__['ir_prob_scaling_factor'] = 'Immune response probability scaling coefficient'
ir_prob_scaling_factor = 1.0 / 100.0
