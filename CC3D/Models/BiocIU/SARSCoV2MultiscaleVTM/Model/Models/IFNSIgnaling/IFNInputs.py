__param_desc__ = {}

__param_desc__['plot_ifn_data_freq'] = 'Plot ifn model data frequency'
plot_ifn_data_freq = 0
__param_desc__['write_ifn_data_freq'] = \
    'Write ifn model data to simulation directory frequency'
write_ifn_data_freq = 0

__param_desc__['plot_med_diff_data_freq'] = 'Plot ifn model medium data frequency'
plot_med_diff_data_freq = 10
__param_desc__['write_med_diff_data_freq'] = \
    'Write ifn model medium data to simulation directory frequency'
write_med_diff_data_freq = 0


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
__param_desc__["k31"] = "IRF7 induction via STATP (1/h)"
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
n = 0.3
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
__param_desc__["IFNe_diffusion_coefficient"] = "IFNe diffusion coefficient (vl^2/s)"
IFNe_diffusion_coefficient = 1.0/10.0
__param_desc__["t2"] = "IFNe_decay (1/hr)"
t2 = 3.481
__param_desc__["Virus_diffusion_coefficient"] = "Virus diffusion coefficient (vl^2/s)"
virus_diffusion_coefficient = 1.0/10.0
__param_desc__["c"] = "IFNe_decay (1/days)"
c = 13.0

# Cell Transition Parameters
__param_desc__['internalization_rate'] = 'Viral Infectivity (PFU/(mL*days))'
b = 2.4E-4
__param_desc__['eclipse_phase'] = '1/days'
k = 4.0