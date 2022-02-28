# coding=utf-8
__param_desc__ = {}

# Data Control
__param_desc__['plot_pop_data_freq'] = 'Plot population data frequency'
plot_pop_data_freq = 0  # Plot population data frequency (disable with 0)
__param_desc__['write_pop_data_freq'] = 'Write population data to simulation directory frequency'
write_pop_data_freq = 1  # Write population data to simulation directory frequency (disable with 0)

__param_desc__['plot_ifn_data_freq'] = 'Plot ifn model data frequency'
plot_ifn_data_freq = 0
__param_desc__['write_ifn_data_freq'] = 'Write ifn model data to simulation directory frequency'
write_ifn_data_freq = 1

__param_desc__['plot_med_diff_data_freq'] = 'Plot ifn model medium data frequency'
plot_med_diff_data_freq = 0
__param_desc__['write_med_diff_data_freq'] = 'Write ifn model medium data to simulation directory frequency'
write_med_diff_data_freq = 1

__param_desc__['plot_plaque_assay_data_freq'] = 'Plot ifn plaque assay data frequency'
plot_plaque_assay_data_freq = 0
__param_desc__['write_plaque_assay_data_freq'] = 'Write ifn plaque assay data to simulation directory frequency'
write_plaque_assay_data_freq = 1

# IFN Model Parameters
__param_desc__["k11"] = "RIGI sensing (μM/h)"
k11 = 0.0
__param_desc__["k12"] = "TLR Sensing (1/h)"
k12 = 9.746
__param_desc__["k13"] = "TLR Sensing (unitless)"
k13 = 12.511
__param_desc__["k14"] = "IFN production via IRF7P (1/h)"
k14 = 13.562
__param_desc__["k21"] = "IFN export rate (1/h)"
k21 = 10.385
__param_desc__["k31"] = "JAK/STAT activation from IFNe (μM/h)"
k31 = 45.922
__param_desc__["k32"] = "JAK/STAT activation from IFNe (μM)"
k32 = 5.464
__param_desc__["k33"] = "JAK/STAT activation from IFNe (unitless)"
k33 = 0.068
__param_desc__["t3"] = "STATP dephosphorylation (1/h)"
t3 = 0.3
__param_desc__["k41"] = "IRF7 induction via STATP (1/h)"
k41 = 0.115
__param_desc__["k42"] = "IRF7 induction via IRF7P (1/h)"
k42 = 1.053
__param_desc__["t4"] = "IRF7 decay (1/h)"
t4 = 0.75
__param_desc__["k51"] = "IRF7 phosphorylation (1/h)"
k51 = 0.202
__param_desc__["t5"] = "IRF7P dephosphorylation (1/h)"
t5 = 0.3
__param_desc__["n"] = "Kinetic order term (unitless)"
n = 3.0
__param_desc__["RIGI"] = "RIGI levels (unitless)"
RIGI = 1.0

# Viral Replication Model Parameters
__param_desc__["k61"] = "Loss of cell viability (1/hr)"
k61 = 0.635
__param_desc__["k71"] = "Viral replication rate (1/μM)"
k71 = 1.537
__param_desc__["k72"] = "Saturation of viral production (1/μM)"
k72 = 47.883
__param_desc__["k73"] = "Viral export rate (1/hr)"
k73 = 0.197

# Diffusible Particles Parameters

__param_desc__['possible_media_for_diffusion'] = "We estimated the diffusion coeficients for these meadia"
possible_media_for_diffusion = ["microcospic mucus", "bulk mucus", "extracellular fluid", "water", "original",
                                "original according to repo"]

__param_desc__['media_selection'] = "Index of media to use"
media_selection = 4  # max 5

__param_desc__["IFNe_diffusion_coefficient"] = "IFNe diffusion coefficient (um^2/min)"
IFNe_diffusion_coefficient = {"microcospic mucus": 0.81 * 60,
                              "bulk mucus": 2.16 * 1e-3 * 60,
                              "extracellular fluid": 4.07 * 60,
                              "water": 46.83 * 60,
                              "original": 3240,  # 54 * 60
                              "original according to repo": 9 / 10}
# as a proportion of the original
# [0.001111111111111111, 2.7962962962962963e-06, 0.005259259259259259, 0.060703703703703704, 1.0, 0.0002777777777777778]

__param_desc__["t2"] = "IFNe decay rate (1/hr)"
t2 = 3.481
__param_desc__["virus_diffusion_coefficient"] = "Virus diffusion coefficient (um^2/min)"
virus_diffusion_coefficient = {"microcospic mucus": 0.06 * 60,
                               "bulk mucus": 1.51 * 1e-4 * 60,
                               "extracellular fluid": 0.284 * 60,
                               "water": 3.278 * 60,
                               "original": 3240,
                               "original according to repo": 9/10}

#

__param_desc__["c"] = "virus decay rate (1/days)"
c = 13.0

# Cell Transition Parameters
__param_desc__['b'] = 'Viral Infectivity (PFU/(mL*days))'
b = 2.4E-4
__param_desc__['k'] = '1/days'
k = 4.0

# death mechanism improvements
__param_desc__["lytic_death_probability"] = "probability of lytic death occuring (i.e., pyroptosis, necroptosis, " \
                                            "simmilar)"

lytic_death_probability = 1

__param_desc__["lytic_death_IFN_release_proportion"] = "proportion of intracellular IFN to be released on lytic death"
lytic_death_IFN_release_proportion = 1


