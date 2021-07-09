import sys
from numpy import log
from .DDMUtils import SetImporter
# import os
# sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
# sys.path.append(os.environ["ViralInfectionVTM"])
# sys.path.insert(1, r'../../Simulation')
# sys.path.append(r'../../Simulation')
# sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from Simulation.ViralInfectionVTMModelInputs import s_to_mcs, um_to_lat_width
except ModuleNotFoundError:
    from ViralInfectionVTMModelInputs import s_to_mcs, um_to_lat_width
# from ViralInfectionVTMModelInputs import s_to_mcs, um_to_lat_width
rate_sets_dict = SetImporter.import_sets_as_dict()

__param_desc__ = {}

# what rates set to use

__param_desc__['set_numb'] = 'chooses which rates to use'
set_numb = 1

set_name = 'set' + str(set_numb)

rs = rate_sets_dict[set_name]  # rates now can be accessed by rs.rate

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
plot_ddm_data_freq = 1  # Plot drug dosing model data frequency (disable with 0)
__param_desc__['write_ddm_data_freq'] = 'Write drug model data to simulation directory frequency'
write_ddm_data_freq = 1  # Write drug dosing model data to simulation directory frequency (disable with 0)

# DDM SBML model

# Treatment options

__param_desc__['intercell_var'] = 'bool for doing inter cell variability of variables'
intercell_var = True

__param_desc__['params_to_var'] = 'which parameters to do the intercell var'
params_to_var = ['rmd_in_rate', 'rmd_out_rate']  # 'protein_rate', 'ace2'

__param_desc__['use_simple_pk'] = 'bool for using simple pk'
use_simple_pk = True

__param_desc__['constant_drug_concentration'] = 'bool flag for constant prodrug'
constant_drug_concentration = False

__param_desc__['prophylactic_treatment'] = 'bool flag for prophylactic treatment'
prophylactic_treatment = False

__param_desc__['treatment_ends'] = 'bool flag for setting a end time for treatment or not'
treatment_ends = False

__param_desc__['sanity_run'] = 'bool for shutting off drug treatment (True -> no treatment)'
sanity_run = False

__param_desc__['double_sbml_step'] = 'bool for doing 2 sbmls calls'
double_sbml_step = False

__param_desc__['double_loading_dose'] = 'bool for having the loading dose be double'
double_loading_dose = True

__param_desc__['first_dose_doubler'] = 'multiplier to double loading dose'
first_dose_doubler = 1
if double_loading_dose:
    first_dose_doubler = 2


__param_desc__['use_alignment'] = 'Bool to start counting time to treatment from N infected cells'
use_alignment = True

__param_desc__['alignt_at_pop'] = 'N infected cells for alignment'
alignt_at_pop = 10

# initial drug concentrations
__param_desc__['Drug_pls'] = 'Concentration of Drug already in plasma'
Drug_pls = 0

__param_desc__['Drug_peri'] = 'Concentration of Drug already in the periphery'
Drug_peri = 0

__param_desc__['Drug_lung'] = 'Concentration of Drug already in the lungs'
Drug_lung = 0

__param_desc__['Ala_met'] = 'Concentration of alanine metabolite already in the system'
Ala_met = 0

__param_desc__['NMP_met'] = 'Concentration of NMP metabolite already in the system'
NMP_met = 0

__param_desc__['NTP_met'] = 'Concentration of NTP metabolite already in the system'
NTP_met = 0

# dosing
__param_desc__['first_dose'] = 'time of first dose in days'
first_dose = 1
if prophylactic_treatment:
    first_dose = 0

__param_desc__['prophylactic_time'] = 'Number of days of prophylactic treatment'
prophylactic_time = 1

if not prophylactic_treatment:
    prophylactic_time = 0

__param_desc__['daily_dose'] = '(twice) how much drug is given in a DAY in mg. Loading dose is double the normal, ' \
                               'it gets galved automatically'
daily_dose = 200

__param_desc__['dose_interval'] = 'time interval between doses in days'
dose_interval = 1

__param_desc__['initial_dose'] = 'initial dose is double the standard [mg]. normal dose is halved in the sbml string'
initial_dose = daily_dose  # / (24. / dose_interval)  # daily_dose/(N doses/day)

__param_desc__['dose'] = 'dose of subsequent treatments [mg]'
dose = initial_dose

if sanity_run:
    initial_dose = 0
    dose = 0

__param_desc__['dose_end'] = 'time of end of treatment in days'
dose_end = 1
if not treatment_ends:
    dose_end = 1e99

__param_desc__['t_half'] = 'half life of the active component (hours)'
t_half = 30.4

__param_desc__['t_half_mult'] = 'half life of the active component (hours)'
t_half_mult = 1

# rate reduction parameters

# parameters

__param_desc__['vary_IC50'] = 'bool to indicate if varying the IC50 or vary a max drug'
vary_IC50 = True

__param_desc__['active_met_ic50'] = 'value for the active metabolite ic50. parameters obtained by running the simple ' \
                                    'PK for 350 hours, detecting the peaks and troths and obtaining the average of ' \
                                    'those [mole/litter]'
active_met_ic50 = 7.89702681970606e-06

__param_desc__['ic50_multiplier'] = 'used to easily change the ic50'
ic50_multiplier = 1

__param_desc__['hill_coef'] = 'Hill coeficient for diminishing hill function of rmax'
hill_coef = 2

__param_desc__['diffusing_drug'] = 'bool to use scalar prodrug or as diffusive species'
diffusing_drug = False

__param_desc__['prodrug_diff_coef_um2'] = 'estimated diffusion length of prodrug (remdesivir) [um^2/s]'
prodrug_diff_coef_um2 = 3

__param_desc__['prodrug_diff_coef_au'] = 'estimated diffusion length of prodrug (remdesivir) [cc3d units]'
prodrug_diff_coef_au = prodrug_diff_coef_um2 * s_to_mcs / (um_to_lat_width ** 2)
