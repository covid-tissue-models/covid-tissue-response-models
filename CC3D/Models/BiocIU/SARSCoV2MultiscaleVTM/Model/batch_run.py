 import os
os.environ["ViralInfectionVTM"] = os.path.dirname(__file__)

# ----------------------------- Setup Instructions ----------------------------- #
# Using this batch script requires three setup steps
# Step 1. CC3D Simulation setup
#   To implement parameter changes automatically in CC3D simulation scripts, do the following
#   Step 1.1. In ViralInfectionVTMSteppables.py, add to the import statements the following line
#       from BatchRun import BatchRunLib
#   Step 1.2. In the __init__ function of the first called steppable, add the following two lines
#       import ViralInfectionVTMModelInputs
#       BatchRunLib.apply_external_multipliers(__name__, ViralInfectionVTMModelInputs)
# Step 2. Environment setup
#   Actually running this script requires an additional step that depends on operating system.
#   You'll modify a script according to your operating system, and then use that script to execute the parameter sweep
#   specified in this script
#       For Windows, use batch_run_win64.bat
#       For OSX, use batch_run_osx.command (not yet implemented)
#       For Linux, use batch_run_linux.sh (not completely tested)
#   Whatever the script, the following step completes setup for execution of this script
#   Step 2.1. In the script that corresponds to your operating system, ensure that Line 4 that sets "PREFIX_CC3D" is
#       the correct path to the root directory of your local CC3D installation
#   If you are using Linux and built from source, then the Python distribution for CC3D may not be at the path specified
#   in batch_run_linux.sh. You can check by verifying that the CC3D installation directory contains the directory
#   "Python37". If it doesn't, locate where Python37 is installed and make the appropriate modifications to
#   batch_run_linux.sh for the variable "PYTHON_INSTALL_PATH".
# You can now execute the script described in Step 2 to run the parameter sweep defined in this script!

# ----------------------------- Script inputs ----------------------------- #
# Multiplier dictionary for parameter sweep to run
#   To simply run multiple replicas of the inputs specified in the simulation scripts, set this to None
#   Each key is a string of a model input variable to multiply in a parameter set of the sweep
#   Each value is a list of multipliers to run
#   Example: to multiply kon by 0.1 and 1, and ir_delay_coeff by 0.1 and 1, do the following
# mult_dict = {'kon': [0.1, 1],
#              'ir_delay_coeff': [0.1, 1]}
#   This will run the following sweep
#     Parameter set 0: {kon *= 0.1, ir_delay_coeff *= 0.1}
#     Parameter set 1: {kon *= 1  , ir_delay_coeff *= 0.1}
#     Parameter set 2: {kon *= 0.1, ir_delay_coeff *= 1}
#     Parameter set 3: {kon *= 1  , ir_delay_coeff *= 1}
mult_dict = None
# mult_dict = {'rel_avail4_EC50': [0.01, 0.1, .25, .5, .75, 1, 1.5, 2, 5, 10]}
# mult_dict = {'rel_avail4_EC50': [1000]}  # for testing

# full_mult_dit = {'first_dose': [0, 6 / 24, 12 / 24, 48 / 24, 72 / 24],  # missing profilaxis. time of 1st dose in days
#                  'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                  # missing continuous dosing. dose interval in days
#                  'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                  'kon': [1 / 4, 1 / 2, 1]}
#
# # NOTE! each batch should be ran twice, as I'm setting numb of rep to 5
# ddm_batch_1 = {'first_dose': [0, 6 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1]}
# ddm_batch_2 = {'first_dose': [12 / 24, 48 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1]}
# ddm_batch_3 = {'first_dose': [72 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1]}
#
# # ______________________________________________
#
# ddm_batch_4 = {'first_dose': [0, 6 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 2]}
# ddm_batch_5 = {'first_dose': [12 / 24, 48 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 2]}
# ddm_batch_6 = {'first_dose': [72 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 2]}
#
# # ______________________________________________
#
# ddm_batch_7 = {'first_dose': [0, 6 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 4]}
# ddm_batch_8 = {'first_dose': [12 / 24, 48 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 4]}
# ddm_batch_9 = {'first_dose': [72 / 24],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 4]}
# # ______________________________________________
# ddm_batch_10 = {'first_dose': [1],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1 / 4]}
# ddm_batch_11 = {'first_dose': [1],
#                'dose_interval': [4 / 24, 6 / 24, 8 / 24, 12 / 24, 1],
#                'rel_avail4_EC50': [0.01, 0.1, .5, .75, 1, 1.25],
#                'kon': [1/2, 1]}
#
# mult_dict = ddm_batch_8

full_cellularized_dict = {}

set_investigation_dict ={'set_numb': [1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                         'drug_ic50': [10, 2, 1, 0.5, 0.1],
                         'daily_dose': [1],
                         'first_dose': [2]}

# Number of replicas to run per parameter set
num_rep = 10
# Number of simulations to run in parallel per parameter set
#   Simulations are implemented in parallel per set of replicas of each parameter set
#   E.g., if running 10 replicas of 2 sets of parameters, then this will run each set of 10 replicas <num_par> at a time
#   You can start by setting this to the number of cores of your machine's CPU

num_par = 1
# Output frequency of model data per simulation replica
# This is passed to the steppables to specify the frequency of writing data generated by the steppables
model_out_freq = 1
# Output frequency of CC3D simulation data per simulation replica
# This is passed to CC3D to specify the frequency of writing spatial data

out_freq = 50
# Sweep output directory
#   If this is set to None, then outputs will be in <root directory>/CallableCoV2VTM
#     Output structure is
#       <sweep_output_folder>/
#           set_0/ (All replica data of parameter set 0 of sweep)
#               ViralInfectionVTMModelParams.csv (input parameter values before variations)
#               Figs/ (Rendered figures of this simulation replica)
#                   Metrics/ (Scalar metrics figures for this simulation replicas)
#                   Spatial/ (Spatial results figures for this simulation replicas)
#               run_0/ (Simulation replica 0 of this parameter set)
#                   CallableSimInputs.csv (variation multipliers used in this simulation replica)
#                   ir_data.csv (raw data of immune response during this simulation replica)
#                   med_diff_data.csv (raw data of diffusion fields during this simulation replica)
#                   pop_data.csv (raw data of population statistics during this simulation replica)
#                   LatticeData/ (raw data of spatial information during this simulation replica)
#               run_1/
#               ...
#           set_1/
#           ...
# sweep_output_folder = os.path.abspath(os.path.join(os.path.splitdrive(os.getcwd())[0], '/DrugDosing_test'))
sweep_output_folder = r'D:\Google Drive IU\phdStuff\covid 19 project\ddm results\new PK\ddm_batch_8'

# Option to execute sweep simulations
#   Set to False to not run simulations
opt_run_sims = False
# Option to render statistics results
#   Set to False to not generate statistics figures
opt_render_stat = True
# Option to render spatial results
#   Set to False to not generate spatial figures
opt_render_spat = False
# Optional dump folder
#   Once all local work is done on a parameter set, the set directory is moved to this location
#   Set to None to leave results where they are first generated
dump_folder = None
# Input modules to export
#   Format input_modules as a list with input modules to export their parameter values as csv
#   Set input_modules to None to not export anything

from Simulation import ViralInfectionVTMModelInputs
from Models.DrugDosingModel import DrugDosingInputs

input_modules = [ViralInfectionVTMModelInputs, DrugDosingInputs]

# ----------------------------- Advanced inputs ----------------------------- #

# Statistics plot manipulators
#   Set to None to pass nothing
#   Manipulators are specified by the items in BatchRun.BatchPostCoV2VTM.export_data_desc
#   Each manipulator is directly called on the returned fig, ax from matplotlib.pyplot.plot()
#   For example, the following sets the vertical limits of the plot for uninfected cells to [0, 900],
#       stat_plot_manips = {'Uninfected': lambda fig, ax: ax.axes.set_ylim(bottom=0, top=900)}
stat_plot_manips = None
stat_plot_manips = {'Uninfected': lambda fig, ax: ax.axes.set_ylim(bottom=0, top=900)}
# ----------------------------- Begin computer work ----------------------------- #
import logging

from nCoVToolkit import nCoVUtils
from BatchRun.CallableCoV2VTM import generic_root_output_folder, CallableCoV2VTMScheduler
from BatchRun.BatchPostCoV2VTM import CoV2VTMSimRunPost, CallableCC3DDataRenderer
from BatchRun.BatchRunLib import cc3d_batch_key


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
            sweep_var_o = sweep_vars[k-1]
            sweep_idx[sweep_var] = int((recur_vals[sweep_var_o] - sweep_idx[sweep_var_o]) / len_mults[sweep_var_o])
        else:
            sweep_var_o = sweep_vars[k-1]
            recur_vals[sweep_var] = int((recur_vals[sweep_var_o] - sweep_idx[sweep_var_o]) / len_mults[sweep_var_o])
            sweep_idx[sweep_var] = recur_vals[sweep_var] % len_mults[sweep_var]

    return {sweep_vars[k]: mult_dict[sweep_vars[k]][sweep_idx[sweep_vars[k]]] for k in range(len(sweep_vars))}


# Example of usage / convenience sequence to do intended overall workflow
if __name__ == '__main__':
    # Check inputs
    if mult_dict is not None:
        assert isinstance(mult_dict, dict), 'mult_dict must be a dictionary or None'

        def check_multipliers(m: dict):
            for var, mults in m.items():
                in_module = False
                for mod in input_modules:
                    if var in mod.__dict__.keys():
                        in_module = True
                        print(f'Found {var} in {mod.__name__}')
                        break
                assert in_module, f'{var} in mult_dict was not found in an input module'
                assert isinstance(mults, list), f'Multipliers for {var} must be specified in a list'
                for x in range(len(mults)):
                    try:
                        mults[x] = float(mults[x])
                    except ValueError:
                        raise ValueError(f'Multiplier {mults[x]} for {var} is not a number')

        check_multipliers(mult_dict)

    assert num_rep > 0, 'Number of replicas per parameter set (num_rep) must be greater than zero'
    if num_par < 1:
        num_par = 1
    assert out_freq > 0, 'Output frequency (out_freq) must be greater than zero'
    if sweep_output_folder is not None:
        assert os.path.isdir(sweep_output_folder), 'Output directory (sweep_output_folder) does not exist'
    if dump_folder is not None:
        assert os.path.isdir(dump_folder), 'Dump directory (dump_folder) does not exist'

    # Setup batch run
    if sweep_output_folder is not None:
        _root_output_folder = sweep_output_folder
    else:
        _root_output_folder = generic_root_output_folder

    num_sets = 1
    if mult_dict is not None:
        for v in mult_dict.values():
            num_sets *= len(v)

    sim_input = list()
    for set_idx in range(num_sets):
        si = sim_input_generator(set_idx)
        si['__param_desc__'] = get_param_descr()
        # Append system configuration inputs
        si[cc3d_batch_key] = {'out_freq': out_freq}
        sim_input.append(si)

    sim_run_sch = CallableCoV2VTMScheduler(root_output_folder=_root_output_folder,
                                           output_frequency=out_freq,
                                           screenshot_output_frequency=out_freq,
                                           num_workers=num_par,
                                           num_runs=num_rep,
                                           sim_input=sim_input,
                                           dump_dir=dump_folder)
    if opt_run_sims:
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

    if opt_render_stat:
        try:
            for set_idx in range(num_sets):
                _cov2_vtm_sim_run = sim_run_sch.run_instance(set_idx)
                if sim_run_sch.is_dumping:
                    _cov2_vtm_sim_run.output_dir_root = sim_run_sch.dump_set_directory(set_idx)
                else:
                    _cov2_vtm_sim_run.output_dir_root = sim_run_sch.output_set_directory(set_idx)
                cov2_vtm_sim_run_post = CoV2VTMSimRunPost(_cov2_vtm_sim_run)
                cov2_vtm_sim_run_post.export_transient_plot_trials(manipulators=stat_plot_manips)
        except Exception as err:
            logging.exception('Error during transient plot rendering.')
            opt_render_stat = False
            print('Disabling transient plot rendering')

    if opt_render_spat:
        data_dirs = []
        out_dirs = []
        set_labs = []
        run_labs = []
        for set_idx in range(num_sets):
            if sim_run_sch.is_dumping:
                set_directory = sim_run_sch.dump_set_directory(set_idx)
                get_run_directory = sim_run_sch.dump_run_directory
            else:
                set_directory = sim_run_sch.output_set_directory(set_idx)
                get_run_directory = sim_run_sch.output_run_directory
            fig_directory = os.path.dirname(set_directory)

            for run_idx in range(sim_run_sch.num_runs[set_idx]):
                data_dirs.append(get_run_directory(set_idx, run_idx))
                out_dirs.append(fig_directory)
                set_labs.append(set_idx)
                run_labs.append(run_idx)

        try:
            callable_cc3d_renderer = CallableCC3DDataRenderer(data_dirs=data_dirs,
                                                              out_dirs=out_dirs,
                                                              set_labs=set_labs,
                                                              run_labs=run_labs,
                                                              num_workers=num_par)
            screenshot_specs = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'screenshots.json')
            callable_cc3d_renderer.load_screenshot_specs(screenshot_specs)
            callable_cc3d_renderer.render_trial_results_par(opts={'log_scale': True,
                                                                  'fixed_caxes': True})
        except Exception as err:
            logging.exception('Error during spatial plot rendering.')
            opt_render_spat = False
            print('Diasbling spatial plot rendering.')
