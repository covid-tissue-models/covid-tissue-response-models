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
#       For Linux, use batch_run_linux.sh (not yet implemented)
#   Whatever the script, the following step completes setup for execution of this script
#   Step 2.1. In the script that corresponds to your operating system, ensure that Line 4 that sets "PREFIX_CC3D" is
#       the correct path to the root directory of your local CC3D installation
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
# Number of replicas to run per parameter set
num_rep = 10
# Number of simulations to run in parallel per parameter set
#   Simulations are implemented in parallel per set of replicas of each parameter set
#   E.g., if running 10 replicas of 2 sets of parameters, then this will run each set of 10 replicas <num_par> at a time
#   You can start by setting this to the number of cores of your machine's CPU
num_par = 1
# Output frequency of simulation data per simulation replica
out_freq = 10
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
sweep_output_folder = None
# Option to render statistics results
#   Set to False to not generate statistics figures
opt_render_stat = True
# Option to render spatial results
#   Set to False to not generate spatial figures
opt_render_spat = True

# ----------------------------- Begin computer work ----------------------------- #
import os
import shutil

from nCoVToolkit import nCoVUtils
from BatchRun.CallableCoV2VTM import simulation_fname, generic_root_output_folder, CoV2VTMSimRun, run_cov2_vtm_sims
from BatchRun.BatchPostCoV2VTM import CallableCC3DRenderer, CoV2VTMSimRunPost
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
            from Simulation import ViralInfectionVTMModelInputs
            for var, mults in m.items():
                assert var in ViralInfectionVTMModelInputs.__dict__.keys(), f'{var} in mult_dict is not a model input'
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

    # Setup batch run
    if sweep_output_folder is not None:
        _root_output_folder = sweep_output_folder
    else:
        _root_output_folder = generic_root_output_folder

    num_sets = 1
    if mult_dict is not None:
        for v in mult_dict.values():
            num_sets *= len(v)

    for set_idx in range(num_sets):
        # Get model inputs
        sim_input = sim_input_generator(set_idx)
        sim_input['__param_desc__'] = get_param_descr()
        # Append system configuration inputs
        sim_input[cc3d_batch_key] = {'out_freq': out_freq}

        _cov2_vtm_sim_run = CoV2VTMSimRun(num_runs=num_rep,
                                          num_workers=num_par,
                                          output_frequency=out_freq,
                                          screenshot_output_frequency=out_freq,
                                          root_output_folder=_root_output_folder,
                                          sim_input=sim_input)

        # Execute batch simulations
        _cov2_vtm_sim_run = run_cov2_vtm_sims(_cov2_vtm_sim_run)

        # Export model parameters
        from Simulation import ViralInfectionVTMModelInputs
        export_file_abs = os.path.join(os.path.abspath(_cov2_vtm_sim_run.output_dir_root),
                                       "ViralInfectionVTMModelParams.csv")
        nCoVUtils.export_parameters(ViralInfectionVTMModelInputs, export_file_abs)

        # Post-process metrics
        if opt_render_stat:
            cov2_vtm_sim_run_post = CoV2VTMSimRunPost(_cov2_vtm_sim_run)
            cov2_vtm_sim_run_post.export_transient_plot_trials()

        # Render field data
        if not opt_render_spat:
            continue

        callable_cc3d_renderer = CallableCC3DRenderer(_cov2_vtm_sim_run)
        screenshot_specs = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'screenshots.json')

        for run_oi in range(_cov2_vtm_sim_run.num_runs):
            callable_cc3d_renderer.load_screenshot_specs(screenshot_specs, run_oi)
            callable_cc3d_renderer.load_trial_results(run_oi)

            # Get results extrema for each field
            min_max_dict = callable_cc3d_renderer.get_results_min_max(run_oi)

            print(min_max_dict)

            # Apply log scale to all field renders
            def gd_manipulator(gd):
                gd.draw_model_2D.clut.SetScaleToLog10()

            callable_cc3d_renderer.load_rendering_manipulator(gd_manipulator, run_oi)

            # Apply fixed legends to all field renders
            def sc_manipulator(scm):
                for field_name, min_max in min_max_dict.items():
                    scm.screenshotDataDict[field_name].metadata['MinRangeFixed'] = True
                    scm.screenshotDataDict[field_name].metadata['MinRange'] = math.ceil(min_max[1] * 10) / 10 / 1E6
                    scm.screenshotDataDict[field_name].metadata['MaxRangeFixed'] = True
                    scm.screenshotDataDict[field_name].metadata['MaxRange'] = math.ceil(min_max[1] * 10) / 10

            callable_cc3d_renderer.load_screenshot_manipulator(sc_manipulator, run_oi)

            # Render
            callable_cc3d_renderer.render_trial_results(run_oi)

        # Move run results to parameter set subdirectory
        _set_dir = os.path.join(_root_output_folder, f'set_{set_idx}')
        if not os.path.isdir(_set_dir):
            os.mkdir(_set_dir)
        cts_src = [x for x in os.listdir(_root_output_folder) if not x.startswith('set_')]
        for ct_src in cts_src:
            shutil.move(os.path.join(_root_output_folder, ct_src), _set_dir)
