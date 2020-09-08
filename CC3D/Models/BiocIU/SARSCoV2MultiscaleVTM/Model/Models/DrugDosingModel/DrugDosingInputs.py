__param_desc__ = {}

# Data control options
__param_desc__['plot_drug_data_freq'] = 'Plot drug model data frequency'
plot_ddm_data_freq = 1  # Plot recovery model data frequency (disable with 0)
__param_desc__['write_drug_data_freq'] = 'Write drug model data to simulation directory frequency'
write_ddm_data_freq = 0  # Write recovery model data to simulation directory frequency (disable with 0)

# parameters

# DDM SBML model
# initial drug concentrations
__param_desc__['Initial_Drug'] = 'Amount of Drug already in the system'
Drug = 0

__param_desc__['Initial_Available1'] = 'Bioavailable drug already in the system'
Available1 = 0

__param_desc__['Initial_Available2'] = 'bioavailable metabolite 2 already in the system'
Available2 = 0

__param_desc__['Initial_Available3'] = 'bioavailable metabolite 3 already in the system'
Available3 = 0

__param_desc__['Initial_Available4'] = 'bioavailable metabolite 4 already in the system'
Available4 = 0

# rates
__param_desc__['k0'] = 'bioavailability rate, units /day. 100/day <-> 15 minutes'
k0 = 100.0

__param_desc__['d0'] = 'clearance time of drug, units /day'
d0 = 1.0

__param_desc__['k1'] = 'metabolism of primary drug rate, units /day'
k1 = 25.0

__param_desc__['d1'] = 'clearance time of avail1, units /day'
d1 = 6.0

__param_desc__['k2'] = 'metabolism of secondary product, units /day'
k2 = 25.0

__param_desc__['d2'] = 'clearance time of avail2, units /day'
d2 = 6.0

__param_desc__['k3'] = 'metabolism of tertiary product, units /day'
k3 = 25.0

__param_desc__['d3'] = 'clearance time of avail3, units /day'
d3 = 6.0

__param_desc__['d4'] = 'clearance time of avail4, units /day'
d4 = 6.0

# dosing
__param_desc__['first_dose'] = 'time of first dose in days'
first_dose = 0.5

__param_desc__['initial_dose'] = 'initial dose (arbitrary amount)'
initial_dose = 10

__param_desc__['dose_interval'] = 'time interval between doses in days'
dose_interval = 0.25

__param_desc__['dose'] = 'dose of subsequent treatments (arbitrary units)'
dose = 10


# rate reduction parameters

__param_desc__['rel_avail4_EC50'] = 'EC50 value for rmax reduction in therms of max available 4,' \
                                    ' ie EC50 = rel_avail4_EC50 * max(available 4)'
rel_avail4_EC50 = 1

__param_desc__['hill_coef'] = 'Hill coeficient for diminishing hill function of rmax'
hill_coef = 2
