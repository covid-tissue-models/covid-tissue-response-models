
# Data control options
track_model_variables = False  # Set to true to enable cell-level tracking of model variables (for rendering)
plot_vrm_data_freq = 1  # Plot viral replication model data frequency (disable with 0)
write_vrm_data_freq = 1  # Write viral replication model data to simulation directory frequency (disable with 0)
plot_vim_data_freq = 0  # Plot viral internalization model data frequency (disable with 0)
write_vim_data_freq = 0  # Write viral internalization model data to simulation directory frequency (disable with 0)
plot_pop_data_freq = 0  # Plot population data frequency (disable with 0)
write_pop_data_freq = 0  # Write population data to simulation directory frequency (disable with 0)
plot_med_diff_data_freq = 0  # Plot total diffusive field amount frequency (disable with 0)
write_med_diff_data_freq = 0  # Write total diffusive field amount frequency (disable with 0)
plot_ir_data_freq = 0  # Plot immune recruitment data frequency (disable with 0)
write_ir_data_freq = 0  # Write immune recruitment data to simulation directory frequency (disable with 0)
plot_spat_data_freq = 0  # Plot spatial data frequency (disable with 0)
write_spat_data_freq = 0  # Write spatial data to simulation directory frequency (disable with 0)

# Conversion Factors
s_to_mcs = 120.0  # s/mcs
um_to_lat_width = 4.0  # um/lattice_length

pmol_to_cc3d_au = 1e15  # 1e15au/1pmol

# Experimental Parameters
exp_cell_diameter = 12.0  # um

exp_unpacking_rate = 1.0 / 10.0 * 1.0 / 60.0
exp_replicating_rate = 1.0 / 20.0 * 1.0 / 60.0  # 1.0/20.0min * 1.0min/60.0s = 1.0/1200.0s
exp_translating_rate = 1.0 / 30.0 * 1.0 / 60.0
exp_packing_rate = 1.0 / 10.0 * 1.0 / 60.0
exp_secretion_rate = 1.0 / 10.0 * 1.0 / 60.0

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

exp_max_ck_diff_len = 100  # um  from (A); this diffusion length is due to uptake, I'm using it as a base
# for the decay due to ck leakage outside of the simulatted lattice.
# dl = sqrt(D/g) -> g = D/dl**2
exp_min_ck_decay = exp_cytokine_dc_cyto/(1.1*exp_max_ck_diff_len)**2

# molecule / (cell second); maximum consumption of cytokine; actually a range [0.3,1] molecule / (cell second) (A)
exp_max_cytokine_consumption = 1
# molecule / (cell second) (B)
exp_max_cytokine_immune_secretion = 10

exp_max_cytokine_consumption_mol = 3.5e-4  # pM/s
exp_max_cytokine_immune_secretion_mol = 3.5e-3  # pM/s

exp_EC50_cytokine_immune = 1  # pM from (B), it's a range from [1,50]pM
# tbd: try to find experimental data
minimum_activated_time_seconds = 60 * 60  # min * s/min

exp_max_amount_viral_mRNA = 1000

# oxidation agent

exp_oxi_dl = 3 * exp_cell_diameter  # [um]; guestimation; [.3,3]
exp_oxi_dc_water = 1.2  # cm2/day; http://www.idc-online.com/technical_references/pdfs/chemical_engineering/Transport_Properties_of_Hydrogen_Peroxide_II.pdf
exp_oxi_dc_water = exp_oxi_dc_water * 1e8 / 86400  # um2/s
exp_oxi_dc_cyto = exp_oxi_dc_water * .16  # rescale by relative density; cyto ~ 6*water
# exp_oxi_dc_cyto = exp_cytokine_dc_cyto

# the experimental values are WAY WAY too high for the simulation to behave properlly, so:
exp_oxi_dc_cyto = 4 * exp_cytokine_dc_cyto


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
ck_equilibrium = 2.1*EC50_ck_immune  # equilibrium amount of ck in immune surface
ck_memory_immune = 1 - max_ck_consume/ck_equilibrium  # decay therm for "seen" ck by immune

max_ck_secrete_infect = 10*max_ck_secrete_im
minimum_activated_time = minimum_activated_time_seconds/s_to_mcs  # mcs

ec50_infecte_ck_prod = 0.1  # amount of 'internal assembled virus' to be at 50% ck production; chosen from
# tipical simulation values of cell.dict['Uptake'] + cell.dict['Assembled']. they stay around .1 and go up as the
# simulation progresses

# oxidation agent


oxi_dl = exp_oxi_dl/um_to_lat_width

oxi_dc = exp_oxi_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)

oxi_decay = oxi_dl**2/oxi_dc

oxi_sec_thr = 10

max_oxi_secrete = max_ck_secrete_infect

oxi_death_thr = 1.5


# Threshold at which cell infection is evaluated
cell_infection_threshold = 1.0

# Threshold at which cell death is evaluated
cell_death_threshold = 1.2
# Probability of survival of infected cell once cell_death_threshold is reached
survival_probability = 0.95
# Hill equations coefficients for probability of viral-induced apoptosis
# Measurements are taken w.r.t. the total amount of assembled viral particles in a cell's simulation subdomain
# dissociationt constant
diss_coeff_uptake_apo = 100.0
# Hill coefficient
hill_coeff_uptake_apo = 2.0

# Hill equation coefficients for probability of viral particle uptake from the environment
# Measurements are taken w.r.t. the total amount of viral particles in a cell's simulation subdomain
# Hill coefficient
hill_coeff_uptake_pr = 2.0

# Efficiency of viral uptake
relative_viral_uptake = 0.1

# Number of immune cells to seed at the beginning of the simulation
# Leave as zero if running through CallableCC3D (bug in NeighborTracker)
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
vi_step_size = vr_step_size

# Viral Internalization parameters
exp_kon = 1.4E5  # 1/(M * s)
exp_koff = 1.4E-4  # 1/s
exp_internalization_rate = 1.0/10.0  # 1/s

# Initial number of receptors on the cell surface
initial_unbound_receptors = 200

# Number of cell receptors at which the replication rate is half max (sets the steady state value of mRNA)
r_half = exp_max_amount_viral_mRNA/(replicating_rate/translating_rate-1)

kon = exp_kon * s_to_mcs * 1.0E15 * (1.0/(um_to_lat_width**3)) * (1.0/1.0E12) * (1.0/pmol_to_cc3d_au) * 100
koff = exp_koff * s_to_mcs
internalization_rate = exp_internalization_rate * s_to_mcs

print('unpacking_rate=', unpacking_rate)
print('replicating_rate=', replicating_rate)
print('translating_rate=', translating_rate)
print('packing_rate=', packing_rate)
print('r_half=', r_half)

# Imunne recruitment parameters

# State variable rate addition coefficient
ir_add_coeff = 1.0
# State variable rate subtraction coefficient
if initial_immune_seeding == 0:
    ir_subtract_coeff = ir_add_coeff / 5.0
else:
    ir_subtract_coeff = ir_add_coeff / initial_immune_seeding
# State variable rate delay coefficient
ir_delay_coeff = 1*1E-3
# State variable rate decay coefficient
ir_decay_coeff = 1E-1
# Ratio of decay cytokine used as recruitment signal
ir_transmission_coeff = 5E-1
# Scales state variable in probability functions
ir_prob_scaling_factor = 1.0 / 100.0
