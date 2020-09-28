from numpy import log

__param_desc__ = {}

# remdesivir metabolites names

__param_desc__['remdesivir_name'] = 'proper molecular name, from https://doi.org/10.1111/cts.12840 '

remdesivir_name = 'RDV;  GS-5734'

__param_desc__['intermediary_metabolite_1'] = 'proper molecular name for intermediary metabolite, ' \
                                              'from https://doi.org/10.1111/cts.12840 '
intermediary_metabolite_1 = 'GS-704277'

__param_desc__['intermediary_metabolite_2'] = 'proper molecular name for intermediary metabolite, ' \
                                              'from https://doi.org/10.1111/cts.12840 '
intermediary_metabolite_2 = 'GS-441524'

__param_desc__['active_met_name'] = 'proper molecular name for active metabolite (tri-phosphate)'
active_met_name = 'GS-443902'

# Data control options
__param_desc__['plot_ddm_data_freq'] = 'Plot drug model data frequency'
plot_ddm_data_freq = 1  # Plot recovery model data frequency (disable with 0)
__param_desc__['write_ddm_data_freq'] = 'Write drug model data to simulation directory frequency'
write_ddm_data_freq = 0  # Write recovery model data to simulation directory frequency (disable with 0)

# parameters
__param_desc__['auto_ec50'] = 'bool for auto scaling of EC50 by max(avail4) and rel_avail4_EC50'
auto_ec50 = False

__param_desc__['ec50'] = 'value for ec50 in the hill equation, only used if auto_ec50 is false'
ec50 = 4.14360796

# DDM SBML model

# Treatment options

constant_drug_concentration = False

profilactic_treatment = False

treatment_ends = False

# initial drug concentrations
__param_desc__['Drug'] = 'Amount of Drug already in the system'
Drug = 0

__param_desc__['Available1'] = 'Bioavailable drug already in the system'
Available1 = 0

__param_desc__['Available2'] = 'bioavailable metabolite 2 already in the system'
Available2 = 0

__param_desc__['Available3'] = 'bioavailable metabolite 3 already in the system'
Available3 = 0

__param_desc__['Available4'] = 'bioavailable metabolite 4 already in the system'
Available4 = 0

# rates
# rates derived from https://doi.org/10.1111/cts.12840
__param_desc__['k0'] = 'bioavailability rate, units /day. 100/day <-> 15 minutes'
k0 = 10.0

__param_desc__['d0'] = 'clearance time of drug, units /day'
d0 = 16.635

__param_desc__['k1'] = 'metabolism of primary drug rate, units /day'
k1 = 1.0

__param_desc__['d1'] = 'clearance time of avail1, units /day'
d1 = 8.317

__param_desc__['k2'] = 'metabolism of secondary product, units /day'
k2 = 989.6

__param_desc__['d2'] = 'clearance time of avail2, units /day'
d2 = 8.317

__param_desc__['k3'] = 'metabolism of tertiary product, units /day'
k3 = 158.4

__param_desc__['d3'] = 'clearance time of avail3, units /day'
d3 = 0.693

__param_desc__['active_metabolite_half_life'] = 'half life of active metabolite, available 4, in days. from ' \
                                                'https://doi.org/10.1111/cts.12840 '
active_metabolite_half_life = 22 / 24  # NOT USED. range from 17.2 to 26.9h

__param_desc__['d4'] = 'clearance time of avail4, units /day'
d4 = 0.693

# dosing
__param_desc__['first_dose'] = 'time of first dose in days'
first_dose = 0.5 / 24
if profilactic_treatment:
    first_dose = 0

__param_desc__['initial_dose'] = 'initial dose (arbitrary amount)'
initial_dose = 10

__param_desc__['dose_interval'] = 'time interval between doses in days'
dose_interval = 0.25

__param_desc__['dose'] = 'dose of subsequent treatments (arbitrary units)'
dose = 10

__param_desc__['dose_end'] = 'time of end of treatment in days'
dose_end = 1
if not treatment_ends:
    dose_end = 1e99

# rate reduction parameters

__param_desc__['rel_avail4_EC50'] = 'EC50 value for rmax reduction in therms of max available 4,' \
                                    ' ie EC50 = rel_avail4_EC50 * max(available 4)'
rel_avail4_EC50 = 1

__param_desc__['hill_coef'] = 'Hill coeficient for diminishing hill function of rmax'
hill_coef = 2
