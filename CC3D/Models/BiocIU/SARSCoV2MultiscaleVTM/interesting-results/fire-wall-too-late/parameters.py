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

# cytokine

cytokine_dc = exp_cytokine_dc_cyto * s_to_mcs / (um_to_lat_width ** 2)  # CK diff cst
# [1/g] = [s] -> [1/g]/s_to_mcs = [s]/[s/mcs] = [mcs]
# -> [g] = [1/s] -> [g]*[s_to_mcs] = [1/s]*[s/mcs] = [1/mcs]
cytokine_field_decay = exp_min_ck_decay * s_to_mcs
# pM = pmol/L = pmol/(10^15 um^3) = 10^-15 pmol/(um^3) = 10^-15 * um_to_lat_width^3 pmol/pixel
# pM/s = pM * s_to_mcs / MCS
max_ck_consume = exp_max_cytokine_consumption_mol * um_to_lat_width ** 3 * s_to_mcs * 1e-15 * pmol_to_cc3d_au  # cc3d_au/(pixel seconds)
max_ck_secrete_im = exp_max_cytokine_immune_secretion_mol * um_to_lat_width ** 3 * s_to_mcs * 1e-15 * pmol_to_cc3d_au  # * cc3d_au/(pixel seconds)
EC50_ck_immune = exp_EC50_cytokine_immune * um_to_lat_width ** 3 * 1e-15 * pmol_to_cc3d_au  # * cc3d_au/pixel
ck_equilibrium = 1.5*EC50_ck_immune # equilibrium amount of ck in immune surface
ck_memory_immune = 1 - max_ck_consume/ck_equilibrium # decay therm for "seen" ck by immune

max_ck_secrete_infect = 10*max_ck_secrete_im


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