# Written by T.J. Sego, Ph.D.

import math
from matplotlib import pyplot as plt
import os
import statistics
import sys
import json


sys.path.append(os.path.dirname(os.path.dirname(__file__)))
import render_defs

sys.path.append(os.path.dirname(__file__))
import render_defs_2_5D_infection_sweep as proj_render_defs

from GridderDealWOtherVars import *

from grid_color_picker_functions import get_sets_of_param_value

# Inputs

dir_spatial_results_root_abs = r'C:\my_results' # This should point to the folder that contains the set_* folders


inputs = get_sim_inputs(dir_spatial_results_root_abs)
rate_mult = 1 / 2
rate = rate_mult * 30.4

first_dose = get_parameters_by_name('first_dose', inputs[0]['__input_dict__'])['first_dose']
    
dir_spatial_results_rel = get_sets_of_param_value('first_dose', first_dose, inputs)

dir_root_output_abs = os.path.join(dir_spatial_results_root_abs, 'figures') # This can point anywhere

num_runs = 8
output_frequency = 1
screenshot_output_frequency = 50
plot_sets = False
plot_grid = True

# the 2 variables that will be the axis of the grid
var1_name = 'dose_interval'
var2_name = 'ic50_multiplier'

x_label_str_var = render_defs.var_labels[var1_name]
y_label_str_var = render_defs.var_labels[var2_name]
rc_params = render_defs.rc_params_multiset.copy()
rc_params['figure.figsize'] = (6.5, 4)
compare_ode = False

grid_plot_pads = {'top': 0.05,
                  'topMinor': 0.045,
                  'bottom': 0.03,
                  'left': 0.045,
                  'right': 0.045,
                  'rightMinor': 0.045}

# Plot manipulators
stat_plot_manips = proj_render_defs.stat_plot_manips.copy()

grid_fig_manips = dict()


def grid_axes_manip(figfig):
    [ax.set_xlabel('') for ax in figfig.axes.flat]
    # https://matplotlib.org/3.1.1/api/_as_gen/matplotlib.pyplot.subplots_adjust.html
    figfig.fig.subplots_adjust(top=0.87, bottom=0.12, hspace=0.4, wspace=0.4)


for k in list(stat_plot_manips.keys()):
    grid_fig_manips[k] = grid_axes_manip

export_fig_height = proj_render_defs.export_fig_height
export_fig_width = proj_render_defs.export_fig_width
params_export = proj_render_defs.params_export

# Load environment
import csv
import shutil

sys.path.append(render_defs.dir_env)
from BatchRun.CallableCoV2VTM import *
from BatchRun.BatchPostCoV2VTM import *


#
# [modify_data_desc(k, v) for k, v in render_defs.export_data_desc_g.items()]
# for k, v in render_defs.fig_save_names_g.items():
#     [modify_fig_save_names(k, kk, vv) for kk, vv in v.items()]
# for k, v in render_defs.y_label_str.items():
#     [modify_y_label_str(k, kk, vv) for kk, vv in v.items()]


#

def get_sim_run_set_num(cov2_vtm_sim_run):
    return int(os.path.basename(cov2_vtm_sim_run.output_dir_root).replace("set_", ""))


def get_error_filename(sim_run_post, data_desc, param_name, err_dir_abs):
    return sim_run_post.generate_transient_plot_trials_filename(
        data_desc, param_name, fig_dir=err_dir_abs).replace("metric_", "error_")


def get_run_variations(cov2_vtm_sim_run):
    from BatchRun import BatchRunLib

    var1 = None
    var2 = None
    with open(os.path.join(cov2_vtm_sim_run.get_run_output_dir(0), r'CallableSimInputs.csv')) as csvf:
        csv_reader = csv.reader(csvf, delimiter=',')
        for row in csv_reader:
            # Catch new format
            if row[0] == BatchRunLib.cc3d_input_key:
                input_dict = eval(row[1])
                if var1_name in input_dict.keys():
                    var1 = int(round(input_dict[var1_name] * 24))
                if var2_name in input_dict.keys():
                    var2 = input_dict[var2_name]
                # if isinstance(var1, float):
                #     var1 /= 10.0
                # if isinstance(var2, float):
                #     var2 /= 10.0
                return var1, var2

            var_name = row[0]
            if var_name == var1_name:
                var1 = float(row[1]) / 10.0
            elif var_name == var2_name:
                var2 = float(row[1]) / 10.0

    return var1, var2


def import_ode_batch_results(_ode_res_dir_abs):
    # Import ODE results
    cov2_vtm_sim_run_o = CoV2VTMSimRun(root_output_folder=_ode_res_dir_abs,
                                       output_frequency=output_frequency,
                                       screenshot_output_frequency=screenshot_output_frequency,
                                       num_runs=1)
    temp_output_dir = cov2_vtm_sim_run_o.get_run_output_dir(0)
    if not os.path.isdir(temp_output_dir):
        os.makedirs(temp_output_dir)
    ode_data_files = []
    ode_data_res = []
    for f in os.listdir(_ode_res_dir_abs):
        if f.endswith('.csv'):
            f_name = f.rstrip('.csv')
            if f_name in export_data_desc.keys():
                ode_data_res.append(f_name)
                ode_data_files.append(f)

    for f in ode_data_files:
        shutil.copyfile(os.path.join(_ode_res_dir_abs, f), os.path.join(temp_output_dir, f))

    _sim_run_post_o = CoV2VTMSimRunPost(cov2_vtm_sim_run_o)

    shutil.rmtree(temp_output_dir)
    return _sim_run_post_o.batch_data_summary


def get_color_vir_thr_start_treat_14_treat(batch_path, set_numb, set_color_file='set_colors.json'):
    ax_color = None
    start = None
    day14 = None
    vir_thr = None
    set_name = f'set_{set_numb}'
    if os.path.isfile(os.path.join(batch_path, set_color_file)):
        with open(os.path.join(batch_path, set_color_file), 'r') as jf:
            d = json.load(jf)
        ax_color = d[set_name][1]
        vir_thr = d[set_name][2]
        start = d[set_name][3]
        day14 = d[set_name][4]
    return ax_color, vir_thr, start, day14


if __name__ == '__main__':
    # Check existence of all rendering output subdirectories
    if not os.path.isdir(dir_root_output_abs):
        os.makedirs(dir_root_output_abs)

    # Grid prep
    batch_data_grid = None
    if plot_grid:
        grid_dir_abs = os.path.join(dir_root_output_abs, "Grids")
        if not os.path.exists(grid_dir_abs):
            os.mkdir(grid_dir_abs)
        batch_data_grid = dict()

    # ODE comparison prep
    sim_run_post_list = []
    var_dict = dict()

    for set_idx in range(len(dir_spatial_results_rel)):

        dir_spatial_result_rel = dir_spatial_results_rel[set_idx]
        # dir_ode_result_rel = dir_ode_results_rel[set_idx]

        results_dir = os.path.join(dir_spatial_results_root_abs, dir_spatial_result_rel)
        # ode_results_dir = os.path.join(os.path.dirname(__file__), dir_ode_result_rel)

        print(f'PostMultiSet processing {results_dir}')

        # batch_data_o = import_ode_batch_results(ode_results_dir)

        fig_dir = os.path.join(dir_root_output_abs, dir_spatial_result_rel)
        if not os.path.isdir(fig_dir):
            os.makedirs(fig_dir)

        # Import spatial results
        cov2_vtm_sim_run_s = CoV2VTMSimRun(root_output_folder=results_dir,
                                           output_frequency=output_frequency,
                                           screenshot_output_frequency=screenshot_output_frequency,
                                           num_runs=num_runs)
        sim_run_post_s = CoV2VTMSimRunPost(cov2_vtm_sim_run_s)
        batch_data_s = sim_run_post_s.batch_data_summary
        var_dict[sim_run_post_s] = get_run_variations(cov2_vtm_sim_run_s)

        # Combine
        for k in batch_data_s.keys():
            # print(k)
            batch_data_s[k].pop("batchMean")
            batch_data_s[k].pop("batchStDev")
            # batch_data_s[k]["ODE"] = batch_data_o[k][0]

        # Trim
        for k in batch_data_s:
            max_step = max([list(batch_data_s[k][t].keys())[-1] for t in batch_data_s[k].keys() if t != "ODE"])
            # batch_data_s[k]["ODE"] = {k: v for k, v in batch_data_s[k]["ODE"].items() if k < max_step}

        # Plot individual sets
        if plot_sets:

            params = dict()
            filenames = dict()
            for k in batch_data_s:
                params[k] = sim_run_post_s.return_param_names(k)
                for p in sim_run_post_s.return_param_names(k):
                    if p in params_export:
                        filenames[p] = sim_run_post_s.generate_transient_plot_trials_filename(k, p, fig_dir=fig_dir)

            render_defs.export_final_plots_trials(batch_data_summary=batch_data_s, params=params, filenames=filenames,
                                                  fig_height=export_fig_height, fig_width=export_fig_width,
                                                  manips=stat_plot_manips, rc_params=render_defs.rc_params_base)

        # Get grid info
        if plot_grid:

            # Get parameter variation for this set
            var_x_val, var_y_val = var_dict[sim_run_post_s]

            if var_x_val is None:
                continue

            # Add set data to grid data
            if var_x_val not in batch_data_grid.keys():
                batch_data_grid[var_x_val] = dict()
            batch_data_grid[var_x_val][var_y_val] = batch_data_s

        if compare_ode:
            sim_run_post_list.append(sim_run_post_s)

    # Build grid
    if plot_grid:
        plt.rcParams.update(rc_params)
        x_vals = list(batch_data_grid.keys())
        x_vals.sort()

        y_vals = list()
        for x in x_vals:
            [y_vals.append(y) for y in batch_data_grid[x] if y not in y_vals]
        y_vals.sort()

        n_cols = len(x_vals)
        n_rows = len(y_vals)

        for d in sim_run_post_s.get_data_descs():
            for p in sim_run_post_s.return_param_names(d):
                print(p)
                if p not in params_export:
                    continue

                fp = render_defs.FigurePack((n_rows, n_cols))

                for j in range(n_cols):
                    for i in range(n_rows):
                        x = x_vals[j]
                        y = y_vals[i]
                        if x not in batch_data_grid.keys() or y not in batch_data_grid[x].keys():
                            continue
                        batch_data_summary = batch_data_grid[x][y]
                        fp.set_subplot_queue((i, j))
                        # generate_transient_plot_trials(batch_data_summary, d, p, fp)
                        # ax = fp.get_subplot((i, j))
                        # ax.lines[-1].set_marker("None")
                        # ax.lines[-1].set_color('k')
                        set_numb = j + i * n_cols

                        ax_color, _, start_time, day14_time = get_color_vir_thr_start_treat_14_treat(
                            dir_spatial_results_root_abs, set_numb)
                        # print(j, i, set_numb, ax_color)
                        render_defs.generate_transient_plot_quantiles(batch_data_summary, d, p, fp, color_ax=ax_color,
                                                                      start_time=start_time, day14_time=day14_time)

                manip = None
                if stat_plot_manips is not None and p in stat_plot_manips.keys():
                    manip = stat_plot_manips[p]

                render_defs.format_grid_fig(figfig=fp, param_vals_x=x_vals, param_vals_y=y_vals,
                                            x_label_str_var=x_label_str_var, y_label_str_var=y_label_str_var,
                                            grid_plot_pads=grid_plot_pads, manip=manip)
                if grid_fig_manips is not None and p in grid_fig_manips.keys():
                    grid_fig_manips[p](fp)

                fig_filename = sim_run_post_s.generate_transient_plot_trials_filename(d, p, fig_dir=grid_dir_abs)
                print('Exporting grid ', fig_filename)
                fp.fig.savefig(fig_filename)

                fig_filename = sim_run_post_s.generate_transient_plot_trials_filename(d, p, fig_dir=grid_dir_abs,
                                                                                      fig_suffix='.eps')
                print('Exporting grid ', fig_filename)
                fp.fig.savefig(fig_filename)

                fig_filename = sim_run_post_s.generate_transient_plot_trials_filename(d, p, fig_dir=grid_dir_abs,
                                                                                      fig_suffix='.pdf')
                print('Exporting grid ', fig_filename)
                fp.fig.savefig(fig_filename)

    # Generate forecast metrics
    if compare_ode:
        err_dict_set = dict()
        err_dict_run = dict()
        err_sum_dict = {sim_run_post: 0 for sim_run_post in sim_run_post_list}
        for sim_run_post in sim_run_post_list:
            for d in sim_run_post.get_data_descs():
                if "ODE" not in sim_run_post.batch_data_summary[d].keys():
                    continue
                run_idx_list = [x for x in sim_run_post.batch_data_summary[d].keys() if x != "ODE"]
                ode_sol = sim_run_post.batch_data_summary[d]["ODE"]
                mcs_max_ode = max(ode_sol.keys())
                for p in sim_run_post.return_param_names(d):
                    nrmse_runs = [0] * len(run_idx_list)
                    for run_idx in run_idx_list:
                        run_res = sim_run_post.batch_data_summary[d][run_idx]
                        mcs_list_s = list(run_res.keys())
                        mcs_list = [x for x in mcs_list_s if x <= min(mcs_max_ode, max(mcs_list_s))]

                        min_run_res = min([run_res[mcs][p] for mcs in mcs_list])
                        max_run_res = max([run_res[mcs][p] for mcs in mcs_list])
                        r = max_run_res - min_run_res
                        if r > 0 and len(mcs_list) > 0:
                            nrmse_runs[run_idx] = math.sqrt(sum([(run_res[mcs][p] - ode_sol[mcs][p]) ** 2.0
                                                                 for mcs in mcs_list]) / len(mcs_list)) / r

                    nrmse_mean = statistics.mean(nrmse_runs)
                    nrmse_stdv = statistics.stdev(nrmse_runs)
                    err_sum_dict[sim_run_post] += sum(nrmse_runs)

                    if d not in err_dict_set.keys():
                        err_dict_set[d] = dict()
                    if p not in err_dict_set[d].keys():
                        err_dict_set[d][p] = dict()
                    err_dict_set[d][p][sim_run_post] = (nrmse_mean, nrmse_stdv)
                    if d not in err_dict_run.keys():
                        err_dict_run[d] = dict()
                    if p not in err_dict_run[d].keys():
                        err_dict_run[d][p] = dict()
                    err_dict_run[d][p][sim_run_post] = {r: nrmse_runs[r] for r in run_idx_list}

        err_dir_abs = os.path.join(dir_root_output_abs, "Errors")

        # Save
        data_set_dict = dict()
        data_run_dict = dict()
        for d in err_dict_run.keys():
            for p in err_dict_run[d].keys():
                for sim_run_post, v in err_dict_run[d][p].items():
                    set_num = get_sim_run_set_num(sim_run_post.cov2_vtm_sim_run)
                    var1, var2 = var_dict[sim_run_post]

                    if set_num not in data_set_dict.keys():
                        data_set_dict[set_num] = {'Total': err_sum_dict[sim_run_post],
                                                  var1_name: var1,
                                                  var2_name: var2}
                    nrmse_mean, nrmse_stdv = err_dict_set[d][p][sim_run_post]
                    data_set_dict[set_num][p + ', M'] = nrmse_mean
                    data_set_dict[set_num][p + ', S'] = nrmse_stdv

                    if set_num not in data_run_dict.keys():
                        data_run_dict[set_num] = dict()
                    for r, vv in v.items():
                        if r not in data_run_dict[set_num].keys():
                            data_run_dict[set_num][r] = {var1_name: var1,
                                                         var2_name: var2}
                        data_run_dict[set_num][r][p] = vv

        for set_num in data_run_dict.keys():
            for r in data_run_dict[set_num].keys():
                s = 0
                for d in err_dict_run.keys():
                    for p in err_dict_run[d].keys():
                        s += data_run_dict[set_num][r][p]
                data_run_dict[set_num][r]['Total'] = s

        fieldnames_set = ['Set']
        [fieldnames_set.append(x) for x in data_set_dict[0].keys()]
        fieldnames_run = ['Set', 'Run']
        [fieldnames_run.append(x) for x in data_run_dict[0][0].keys()]
        import csv

        err_rep_filename = os.path.join(os.path.dirname(__file__), f'error_report.csv')
        print('Writing error summary ', err_rep_filename)
        with open(err_rep_filename, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames_set)
            writer.writeheader()
            for k, v in data_set_dict.items():
                row_dict = {kk: vv for kk, vv in data_set_dict[k].items()}
                row_dict['Set'] = k
                writer.writerow(row_dict)

        with open(err_rep_filename, 'a', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames_run)
            writer.writeheader()
            for s in data_run_dict.keys():
                for r in data_run_dict[s].keys():
                    row_dict = {p: v for p, v in data_run_dict[s][r].items()}
                    row_dict['Set'] = s
                    row_dict['Run'] = r
                    writer.writerow(row_dict)

        # Plot
        if not os.path.isdir(err_dir_abs):
            os.makedirs(err_dir_abs)
        for d in err_dict_set.keys():
            for p in err_dict_set[d].keys():
                err_p = {}
                fig_filename = None
                for sim_run_post, v in err_dict_set[d][p].items():
                    nrmse_mean, nrmse_stdv = v
                    set_num = get_sim_run_set_num(sim_run_post.cov2_vtm_sim_run)
                    err_p[set_num] = v
                    fig_filename = get_error_filename(sim_run_post, d, p, err_dir_abs)

                if fig_filename is not None:
                    print('Generating error plot', fig_filename)

                    x_vals = []
                    y_vals = []
                    yerr_vals = []
                    for k, v in err_p.items():
                        x_vals.append(f"Set {k}")
                        y_vals.append(v[0])
                        yerr_vals.append(v[1])

                    f = plt.figure()
                    plt.bar(x_vals, y_vals, yerr=yerr_vals)
                    plt.grid(b='on')
                    plt.xticks(rotation=90)
                    plt.ylabel('Mean RSME')
                    plt.title(y_label_str[d][p])
                    plt.tight_layout()
                    plt.savefig(fig_filename)

        fig_filename = os.path.join(err_dir_abs, "error_total.png")
        print('Generating error plot', fig_filename)
        x_vals = []
        y_vals = []
        for sim_run_post, v in err_sum_dict.items():
            set_num = int(os.path.basename(sim_run_post.cov2_vtm_sim_run.output_dir_root).replace("set_", ""))
            x_vals.append(f"Set {set_num}")
            y_vals.append(v)

        f = plt.figure()
        plt.bar(x_vals, y_vals)
        plt.grid(b='on')
        plt.xticks(rotation=90)
        plt.ylabel('Mean RSME')
        plt.title('Total')
        plt.tight_layout()
        plt.savefig(fig_filename)
