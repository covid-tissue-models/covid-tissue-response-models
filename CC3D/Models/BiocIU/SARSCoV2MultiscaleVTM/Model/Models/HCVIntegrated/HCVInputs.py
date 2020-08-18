__param_desc__ = {}

__param_desc__['plot_ihcv_data_freq'] = 'Plot integrated HCV viral replication model data frequency'
plot_ihcv_data_freq = 0
__param_desc__['write_ihcv_data_freq'] = \
    'Write integrated HCV viral replication model data to simulation directory frequency'
write_ihcv_data_freq = 0

__param_desc__["k1"] = "Tc formation (1/h)"
k1 = 80
__param_desc__["k2"] = "Nascent NS polyprotein translation (1/h)"
k2 = 100
__param_desc__["kc"] = "Viral polyprotein cleavage (1/h)"
kc = 0.6
__param_desc__["kpin"] = "R transport into VMS (1/h)"
kpin = 0.2
__param_desc__["kpout"] = "RP transport into cytoplasm (1/h)"
kpout = 0.2
__param_desc__["kein"] = "ECYT transport into VMS (1/h)"
kein = 1.3E-5
__param_desc__["k3"] = "RIP formation (1/h)"
k3 = 0.02
__param_desc__["k4p"] = "RP synthesis (1/h)"
k4p = 1.7
__param_desc__["k4m"] = "RDS synthesis (1/h)"
k4m = 1.7
__param_desc__["k5"] = "RIDS formation (1/h)"
k5 = 4.0
__param_desc__["upcyt"] = "R degradation (1/h)"
upcyt = 10
__param_desc__["up"] = "RP degradation (1/h)"
up = 0.07
__param_desc__["uds"] = "RDS degradation (1/h)"
uds = 0.06
__param_desc__["uip"] = "RIP degradation (1/h)"
uip = 0.04
__param_desc__["uids"] = "RIDS degradation (1/h)"
uids = 0.13
__param_desc__["utc"] = "TC degradation (1/h)"
utc = 0.015
__param_desc__["ue"] = "E degradation (1/h)"
ue = 0.04
__param_desc__["uecyt"] = "ECYT degradation (1/h)"
uecyt = 0.06

# Initial conditions
__param_desc__["RIBOTOT"] = ""
RIBOTOT = 700.0
__param_desc__["init_rna"] = "Initial number of RNA molecules in initially infected cells"
init_rna = 500

# Conversion factor
__param_desc__["virus_from_ul"] = "Conversion from unitless original model quantities to HCV model units"
virus_from_ul = 10
