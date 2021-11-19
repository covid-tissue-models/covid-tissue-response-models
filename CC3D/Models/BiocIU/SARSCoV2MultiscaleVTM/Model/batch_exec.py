# Model inputs

from investigation_dictionaries import *

# this script takes care of generating jobs and scheduling them using slurm.

# The dictionaries used are defined in investigation_dictionaries.py and are imported here.
# If you want to, e.g., simulate treatment starting with the infection of 10 cells with the half life
# of the antiviral halved you'd do

mult_dict = treatment_starts_0_halved_half_life

num_rep = 5
# Model output frequency
model_out_freq = 1
# Output frequency of simulation data per simulation replica
out_freq = 250

# Root output directory, change to some directory that makes sense to you
sweep_output_folder = r'/home/my_user/my_files'


# Input modules
from Simulation import ViralInfectionVTMModelInputs
from Models.DrugDosingModel import DrugDosingInputs

input_modules = [ViralInfectionVTMModelInputs, DrugDosingInputs]
# Automatic inputs
from BatchRun import BatchRunLib

BatchRunLib.register_auto_inputs(input_module_name='ViralInfectionVTMModelInputs',
                                 plot_var_names=['plot_vrm_data_freq', 'plot_vrm_data_freq', 'plot_vim_data_freq',
                                                 'plot_pop_data_freq', 'plot_ir_data_freq', 'plot_med_diff_data_freq',
                                                 'plot_spat_data_freq', 'plot_death_data_freq'],
                                 write_var_names=['write_pop_data_freq', 'write_med_diff_data_freq',
                                                  'write_ir_data_freq', 'write_death_data_freq'])
BatchRunLib.register_auto_inputs(input_module_name='Models.DrugDosingModel.DrugDosingInputs',
                                 plot_var_names=['plot_ddm_data_freq'],
                                 write_var_names=['write_ddm_data_freq'])

# Carbonate configuration
from BatchRun.BatchRunPrototyping import carbonate_config_template
carbonate_config_template = carbonate_config_template()
carbonate_config_template['jn'] = 'ddm_new_pk_set_inv1'
carbonate_config_template['wh'] = 24
carbonate_config_template['wm'] = 0
carbonate_config_template['ppn'] = 1
carbonate_config_template['mem'] = 10
carbonate_config_template['p'] = 'general'

# Begin computing work
import os
from nCoVToolkit import nCoVUtils
from BatchRun import BatchRunPrototyping

BatchRunPrototyping.simulation_fname = os.path.join(os.path.dirname(__file__), 'ViralInfectionVTM.cc3d')

from BatchRun.BatchRunLib import cc3d_batch_key


def get_num_sets():
    num_sets = 1
    if mult_dict is not None:
        for v in mult_dict.values():
            num_sets *= len(v)
    return num_sets


def get_param_descr():
    if isinstance(mult_dict, dict):
        return {k + "_multiplier": f"Multiplier for {k}" for k in mult_dict.keys()}
    else:
        return None


def sim_input_generator(_set_idx):
    if mult_dict is None:
        return dict()

    sweep_vars = list(mult_dict.keys())
    sweep_idx = {s: 0 for s in sweep_vars}
    len_mults = {s: len(mult_dict[s]) for s in sweep_vars}
    recur_vals = {s: 0 for s in sweep_vars}
    for k in range(len(sweep_vars)):
        sweep_var = sweep_vars[k]
        if k == 0:
            recur_vals[sweep_var] = _set_idx
            sweep_idx[sweep_var] = _set_idx % len_mults[sweep_var]
        elif k == len(sweep_vars) - 1:
            sweep_var_o = sweep_vars[k - 1]
            sweep_idx[sweep_var] = int((recur_vals[sweep_var_o] - sweep_idx[sweep_var_o]) / len_mults[sweep_var_o])
        else:
            sweep_var_o = sweep_vars[k - 1]
            recur_vals[sweep_var] = int((recur_vals[sweep_var_o] - sweep_idx[sweep_var_o]) / len_mults[sweep_var_o])
            sweep_idx[sweep_var] = recur_vals[sweep_var] % len_mults[sweep_var]

    return {sweep_vars[k]: mult_dict[sweep_vars[k]][sweep_idx[sweep_vars[k]]] for k in range(len(sweep_vars))}


def main():
    num_sets = get_num_sets()

    # Add set labels to job names; scheduler adds run labels
    carbonate_config = [carbonate_config_template.copy() for _ in range(num_sets)]
    set_idx = 0
    for cc in carbonate_config:
        cc['jn'] += f'_{set_idx}'
        set_idx += 1

    sim_input = list()
    from BatchRun.CallableCoV2VTM import cc3d_input_key
    for set_idx in range(num_sets):
        si = sim_input_generator(set_idx)
        si['__param_desc__'] = get_param_descr()
        # Append system configuration inputs
        si[cc3d_batch_key] = {'out_freq': model_out_freq}
        sim_input.append({cc3d_input_key: si})

    from BatchRun.BatchRunPrototyping import CallableCC3DCarbonateScheduler
    sim_run_sch = CallableCC3DCarbonateScheduler(carbonate_config=carbonate_config,
                                                 root_output_folder=sweep_output_folder,
                                                 output_frequency=out_freq,
                                                 screenshot_output_frequency=out_freq,
                                                 num_runs=num_rep,
                                                 sim_input=sim_input)
    sim_run_sch.run()

    # Export model parameters
    if input_modules is not None and isinstance(input_modules, list):
        for set_idx in range(num_sets):
            if sim_run_sch.is_dumping:
                o = sim_run_sch.dump_set_directory(set_idx)
            else:
                o = sim_run_sch.output_set_directory(set_idx)

            for x in input_modules:
                export_file_rel = x.__name__.split('.')[-1] + "Params.csv"
                export_file_abs = os.path.join(o, export_file_rel)
                nCoVUtils.export_parameters(x, export_file_abs)


if __name__ == '__main__':
    main()
