# Written by T.J. Sego, Ph.D.

from matplotlib import pyplot as plt
from matplotlib import ticker
import os
from numpy import quantile
from os.path import dirname

# dir_env = os.path.join(dirname(dirname(__file__)), os.path.join('Model', 'GlazierModel', 'Model'))
dir_env = dirname(__file__)
dir_res = os.path.dirname(__file__)

# Plot labels
first_dose_str = 'Time of first dose (days)'
dose_interval_str = 'First dose 1 day after the infection of 10 cells\nTime between doses (hours)'
# dose_interval_str = 'First dose with the infection of 10 cells\nTime between doses (hours)'
rel_avail4_EC50_str = r'Relative IC50 $\kappa$'
kon_str = 'Virus receptor affinity'
ic50_multiplier_str = 'IC50 multiplier'

var_labels = {'first_dose': first_dose_str,
              'dose_interval': dose_interval_str,
              'rel_avail4_EC50': rel_avail4_EC50_str,
              'kon': kon_str,
              'ic50_multiplier': ic50_multiplier_str}

# rcParams definitions
rc_params_base = {'font.family': 'arial',
                  'savefig.dpi': 600}

rc_params_multiset = {'axes.labelsize': 6,
                      'axes.titlesize': 8,
                      'lines.markersize': 2}
rc_params_multiset.update(rc_params_base)

# Data set defs
export_data_desc = {'ir_data': ['ImmuneResp'],
                    'med_diff_data': ['MedViral',
                                      'MedCyt',
                                      'MedOxi'],
                    'pop_data': ['Uninfected',
                                 'Infected',
                                 'InfectedSecreting',
                                 'Dying',
                                 'ImmuneCell',
                                 'ImmuneCellActivated'],
                    'spat_data': ['DeathComp',
                                  'InfectDist'],
                    'death_data': ['Viral',
                                   'OxiField',
                                   'Contact',
                                   'Bystander'],
                    'ddm_rmax_data': ['r_max'],
                    # 'ddm_data': ['Prodrug',
                    #              'ActiveMetabolite'],
                    'ddm_tot_RNA_data': ['Total_viral_RNA_in_cells'],
                    'ddm_mean_RNA_data': ['Mean_viral_RNA_in_cells'],
                    'ddm_total_viral_production_data': ['vir_auc']
                    }

# Results labels
y_label_str = {'ir_data': {'ImmuneResp': 'Immune response state variable'},
               'med_diff_data': {'MedViral': 'Total diffusive virus ($log_{10}$)',
                                 'MedCyt': 'Total diffusive cytokine ($log_{10}$)',
                                 'MedOxi': 'Total oxidative agent'},
               'pop_data': {'Uninfected': 'Number of uninfected cells',
                            'Infected': 'Number of infected cells',
                            'InfectedSecreting': 'Number of infected secreting cells',
                            'Dying': 'Number of dead cells',
                            'ImmuneCell': 'Number of immune cells',
                            'ImmuneCellActivated': 'Number of activated immune cells'},
               'spat_data': {'DeathComp': 'Cell death compactness (ul)',
                             'InfectDist': 'Infection distance (px)'},
               'death_data': {'Viral': 'Number of virally-induced apoptosis deaths',
                              'OxiField': 'Number of oxidative deaths',
                              'Contact': 'Number of cytotoxic kill deaths',
                              'Bystander': 'Number of bystander effect deaths'},
               'ddm_rmax_data': {'r_max': r'$r_{max}$ value'},
               'ddm_data': {'Drug': 'Drug concentration ($log_{10}$)',
                            'ActiveMetabolite': 'Metabolite concentration ($log_{10}$)'},
               'ddm_tot_RNA_data': {'Total_viral_RNA_in_cells': 'Total viral RNA in cells ($log_{10}$)'},
               'ddm_mean_RNA_data': {'Mean_viral_RNA_in_cells': 'Mean viral RNA in cells ($log_{10}$)'},
               'ddm_total_viral_production_data': {'vir_auc': 'Total production of diffusive virus ($log_{10}$)'}
               }

# Figure file names
fig_save_names_g = {'ir_data': {'ImmuneResp': 'metric_immune_response_svar'},
                    'med_diff_data': {'MedViral': 'metric_diffusive_virus',
                                      'MedCyt': 'metric_diffusive_cytokine',
                                      'MedOxi': 'metric_diffusive_oxidator'},
                    'pop_data': {'Uninfected': 'metric_num_uninfected',
                                 'Infected': 'metric_num_infected',
                                 'InfectedSecreting': 'metric_num_infectedSecreting',
                                 'Dying': 'metric_num_dying',
                                 'ImmuneCell': 'metric_num_immune',
                                 'ImmuneCellActivated': 'metric_num_immuneActivated'},
                    'spat_data': {'DeathComp': 'metric_death_compact',
                                  'InfectDist': 'metric_infection_distance'},
                    'death_data': {'Viral': 'metric_death_viral',
                                   'OxiField': 'metric_death_oxi',
                                   'Contact': 'metric_death_contact',
                                   'Bystander': 'metric_death_bystander'},
                    'ddm_rmax_data': {'r_max': 'metric_rmax'},
                    # 'ddm_data': {'Prodrug': 'metric_drug',
                    #              'ActiveMetabolite': 'metric_metabolite'},
                    'ddm_tot_RNA_data': {'Total_viral_RNA_in_cells': 'metric_total_RNA'},
                    'ddm_mean_RNA_data': {'Mean_viral_RNA_in_cells': 'metric_mean_RNA'},
                    'ddm_total_viral_production_data': {'vir_auc': 'metric_vir_AUC'}
                    }

# Plot manipulators
time_label = 'Simulation time (days)'


def manip_plot_log(fig, ax):
    ax.set_yscale('log')


def manip_plot_axis_conv(fig, ax):
    # Convert time units to minutes and relabel
    _scale = 5 * 60 / 60 / 60 / 24  # 5 minutes per step
    ticks = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x * _scale))
    ax.xaxis.set_major_formatter(ticks)


def manip_plot_time_label(fig, ax):
    ax.set_xlabel(time_label)


def manip_set_ticks(xticks=None, yticks=None):
    fx, fy = None, None
    if xticks is not None:
        def fx(fig, ax):
            ax.set_xticks(xticks)
    if yticks is not None:
        def fy(fig, ax):
            ax.set_yticks(yticks)

    def f(fig, ax):
        if fx is not None:
            fx(fig, ax)
        if fy is not None:
            fy(fig, ax)

    return f


def manip_plot_ode(fig, ax):
    ax.lines[-1].set_marker("None")
    ax.lines[-1].set_color('k')


def plot_lim_manip(bottom, top):
    def f(fig, ax):
        ax.axes.set_ylim(bottom=bottom, top=top)

    return f


def plot_hlim_manip(left, right):
    def f(fig, ax):
        ax.axes.set_xlim(left=left, right=right)

    return f


# Prototype addenda to batch post-processing

#
# def generate_transient_subplot_trials(batch_data_summary, data_desc, var_name, fig, ax, manip=None):
#     from BatchRun.BatchPostCoV2VTM import x_label_str_transient
#     # Like generate_transient_plot_trials, but to a passed axis
#     ax.grid()
#
#     data_dict = batch_data_summary[data_desc]
#     for trial_idx in data_dict.keys():
#         if trial_idx in ['batchMean', 'batchStDev']:
#             continue
#         sim_mcs = list(data_dict[trial_idx].keys())
#         y_data = [data_dict[trial_idx][this_mcs][var_name] for this_mcs in sim_mcs]
#         ax.plot(sim_mcs, y_data, label='Trial {}'.format(trial_idx), marker='.')
#
#     ax.set_xlabel(x_label_str_transient)
#     ax.set_ylabel(y_label_str[data_desc][var_name])
#
#     if manip is not None:
#         manip(fig, ax)


color_blue = '#2c7bb6'
color_lblue = '#abd9e9'
color_ltan = '#ffffbf'
color_tan = '#fdae61'
color_red = '#d7191c'

color_blue = '#2c7bb6'
color_lblue = '#abd9e9'
color_ltan = '#ffffbf'
color_tan = '#fdae61'
color_red = '#d7191c'


def generate_transient_subplot_quantiles(batch_data_summary, data_desc, var_name, fig, ax, manip=None):
    from BatchRun.BatchPostCoV2VTM import x_label_str_transient
    # Like generate_transient_plot_trials, but to a passed axis
    ax.grid()

    data_dict = batch_data_summary[data_desc]
    quantiles = [0, 10, 25, 50, 75, 90, 100]
    data_analysis = {}
    data_analyzed = {q: {} for q in quantiles}
    for trial_idx in data_dict.keys():
        if trial_idx in ['batchMean', 'batchStDev', 'ODE']:
            continue
        for m in data_dict[trial_idx].keys():
            if m not in data_analysis.keys():
                data_analysis[m] = []
            data_analysis[m].append(data_dict[trial_idx][m][var_name])

    # ode_dict = None
    # if 'ODE' in data_dict.keys():
    #     ode_dict = data_dict['ODE']

    for q in quantiles:
        for mcs, data_list in data_analysis.items():
            data_analyzed[q][mcs] = quantile(data_list, q / 100)

    roil = data_analyzed[0]
    roiu = data_analyzed[10]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_blue)
    roil = data_analyzed[90]
    roiu = data_analyzed[100]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_blue)

    roil = data_analyzed[10]
    roiu = data_analyzed[25]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_tan)
    roil = data_analyzed[75]
    roiu = data_analyzed[90]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_tan)

    roil = data_analyzed[25]
    roiu = data_analyzed[50]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_lblue)
    roil = data_analyzed[50]
    roiu = data_analyzed[75]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_lblue)

    roi = data_analyzed[50]
    ax.plot(list(roi.keys()), list(roi.values()), color='0', linewidth=0, markersize=1, marker='.')

    # if ode_dict is not None:
    #     roi = {k: ode_dict[k][var_name] for k in ode_dict.keys()}
    #     ax.plot(list(roi.keys()), list(roi.values()), color=color_red, markersize=1, linewidth=1, linestyle='--')

    ax.set_xlabel(x_label_str_transient)
    ax.set_ylabel(y_label_str[data_desc][var_name])

    if manip is not None:
        manip(fig, ax)


generate_transient_subplot_trials = generate_transient_subplot_quantiles


def generate_transient_plot_quantiles(batch_data_summary, data_desc, var_name, fig_pack=None, color_ax=None,
                                      vir_thr=None, start_time=None, day14_time=None):
    from BatchRun.BatchPostCoV2VTM import x_label_str_transient

    fig = None
    if fig_pack is None:
        fig = plt.figure()
        ax = fig.add_subplot(111)
    else:
        ax = fig_pack.get_subplot()


    # if vir_thr is not None:
    #     xmin, xmax = ax.get_xlim()
    #     ax.hlines(vir_thr, xmin=xmin, xmax=xmax, colors='red', linestyles='dotted')

    data_dict = batch_data_summary[data_desc]
    # ode_dict = None
    # if 'ODE' in data_dict.keys():
    #     ode_dict = data_dict['ODE']
    trials_dict = {k: v for k, v in data_dict.items() if k not in ['batchMean', 'batchStDev', 'ODE']}

    quantiles = [0, 10, 25, 50, 75, 90, 100]
    data_analysis = {}
    data_analyzed = {q: {} for q in quantiles}
    for trial_idx in trials_dict.keys():
        for m in data_dict[trial_idx].keys():
            if m not in data_analysis.keys():
                data_analysis[m] = []
            data_analysis[m].append(trials_dict[trial_idx][m][var_name])

    for q in quantiles:
        for mcs, data_list in data_analysis.items():
            data_analyzed[q][mcs] = quantile(data_list, q / 100)

    roil = data_analyzed[0]
    roiu = data_analyzed[10]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_blue)
    roil = data_analyzed[90]
    roiu = data_analyzed[100]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_blue)

    roil = data_analyzed[10]
    roiu = data_analyzed[25]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_tan)
    roil = data_analyzed[75]
    roiu = data_analyzed[90]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_tan)

    roil = data_analyzed[25]
    roiu = data_analyzed[50]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_lblue)
    roil = data_analyzed[50]
    roiu = data_analyzed[75]
    ax.fill_between(list(roil.keys()), list(roil.values()), list(roiu.values()), alpha=0.5, color=color_lblue)

    roi = data_analyzed[50]
    ax.plot(list(roi.keys()), list(roi.values()), color='0', linewidth=0, markersize=1, marker='.')

    # print(start_time)

    # # print(day14_time)
    # if day14_time is not None:
    #     ymin, ymax = ax.get_xlim()
    #     ax.vlines(day14_time, ymin=ymin, ymax=ymax, colors='red', linestyles='solid')

    # if ode_dict is not None:
    #     roi = {k: ode_dict[k][var_name] for k in ode_dict.keys()}
    #     ax.plot(list(roi.keys()), list(roi.values()), color=color_red, markersize=1, linewidth=1, linestyle='--')

    ax.set_xlabel(x_label_str_transient)
    ax.set_ylabel(y_label_str[data_desc][var_name])

    if color_ax is not None:
        ax.grid(True, color=color_ax, alpha=0.5)
        # ax.tick_params(axis='x', colors=color_ax, which='both')
        # ax.tick_params(axis='y', colors=color_ax, which='both')
        for tick in ax.get_yticklines():
            tick.set_color(color_ax)
        for tick in ax.get_xticklines():
            tick.set_color(color_ax)
        for minortick in ax.xaxis.get_minorticklines():
            minortick.set_color(color_ax)
        for side in ['top', 'bottom', 'left', 'right']:
            ax.spines[side].set_color(color_ax)
            ax.spines[side].set_linewidth(1.05)
    else:
        ax.grid()

    # if start_time is not None:
    #     ymin, ymax = ax.get_ylim()
    #     ax.vlines(start_time, ymin=ymin, ymax=ymax, colors='black', linestyles='solid')
    if fig_pack is None:
        fig.tight_layout()
        return fig, ax


def export_final_plots_trials(batch_data_summary, params: dict, filenames: dict, fig_height, fig_width,
                              manips: dict = None, rc_params: dict = None, big_fig_filename: str = None):
    if manips is not None:
        for p in manips.keys():
            assert any([p in pv for pv in params.values()]), f'Manipulator for {p} has not data'

    for k in batch_data_summary:
        assert k in params.keys(), f'No matching parameter list for data set {k}'

    for f in filenames.values():
        assert os.path.isdir(os.path.dirname(f)), f'Containing directory not found for {f}'

    if rc_params is not None:
        plt.rcParams.update(rc_params)

    n_plots = len(filenames.keys())

    big_fig, big_axes = plt.subplots(nrows=1, ncols=n_plots)
    big_fig.set_size_inches(w=n_plots * fig_width, h=fig_height)

    col_idx = 0
    saved_params = list()
    for k in batch_data_summary:
        for p in params[k]:
            if p not in filenames.keys():
                continue
            manip = None
            if manips is not None and p in manips.keys():
                manip = manips[p]
            generate_transient_subplot_trials(batch_data_summary, k, p, big_fig, big_axes[col_idx], manip)
            col_idx += 1
            saved_params.append(p)

    big_fig.tight_layout()
    if big_fig_filename is None:
        big_fig_path = os.path.join(os.path.dirname(__file__), "__bigfig.png")
        keep_big_fig = False
    else:
        big_fig_path = big_fig_filename
        keep_big_fig = True
    print(big_fig_path)
    big_fig.savefig(big_fig_path, bbox_inches='tight')
    plt.close()

    im = plt.imread(big_fig_path)
    h = int(im.shape[0])
    dw = int(im.shape[1] / n_plots)

    fig_counter = 0
    for iw in range(n_plots):
        imc = im[int(0):h, int(iw * dw):int((iw + 1) * dw)]
        fig = plt.figure(figsize=(fig_width, fig_height), dpi=600)
        fig.set_size_inches(w=fig_width, h=fig_height)
        ax = fig.add_subplot(111)
        ax.imshow(imc)
        plt.axis('off')
        fig_filename = filenames[saved_params[fig_counter]]
        print(fig_filename)
        fig.savefig(fig_filename)
        plt.close(fig)
        fig_counter += 1

    if not keep_big_fig:
        os.remove(big_fig_path)


class FigurePack:
    def __init__(self, subplot_specs=None):
        if subplot_specs is None:
            subplot_specs = [0, 0]
        self.fig, self.axes = plt.subplots(subplot_specs[0], subplot_specs[1])

        self.subplot_queue = None
        self.n_rows, self.n_cols = 0, 0
        if subplot_specs is not None:
            self.n_rows = subplot_specs[0]
            self.n_cols = subplot_specs[1]

    def figure(self):
        return self.fig

    def set_subplot_queue(self, _queue_specs):
        assert len(_queue_specs) == 2
        self.subplot_queue = _queue_specs
        if _queue_specs[0] > self.n_rows:
            self.n_rows = _queue_specs[0]
        if _queue_specs[1] > self.n_cols:
            self.n_cols = _queue_specs[1]

    def get_subplot(self, _queue_specs=None):
        if not hasattr(self.axes, 'shape'):
            return self.axes
        if _queue_specs is None:
            _queue_specs = self.subplot_queue
        if (self.axes.shape) == 1:
            return self.axes[max(_queue_specs[0], _queue_specs[1])]
        else:
            return self.axes[_queue_specs[0], _queue_specs[1]]

    def kill_figure(self):
        plt.close(self.fig)
        self.fig = None

    def show(self):
        assert self.fig is not None
        self.fig.show()


def format_grid_fig(figfig: FigurePack, param_vals_x, param_vals_y,
                    x_lims=None, y_lims=None, x_ticks=None, y_ticks=None,
                    x_label_str_var=None, y_label_str_var=None,
                    grid_plot_pads=None, manip=None):
    n_rows = figfig.n_rows
    n_cols = figfig.n_cols

    if manip is not None:
        [manip(figfig.fig, ax) for ax in figfig.axes.flat]

    # Share vertical
    axL = None
    for i in range(n_rows):
        for j in range(n_cols):
            ax = figfig.get_subplot((i, j))

            if j == 0:
                axL = ax
            else:
                axL.get_shared_y_axes().join(axL, ax)

    # Share horizontal axes
    axB = None
    for j in range(n_cols):
        for i in range(n_rows - 1, -1, -1):
            ax = figfig.get_subplot((i, j))

            if i == n_rows - 1:
                axB = ax
            else:
                axB.get_shared_x_axes().join(axB, ax)

    # Set labeling to only outer subplots
    for ax in figfig.axes.flat:
        ax.label_outer()

    # Set common x label
    if n_rows > 1:
        # Set common x label
        x_lab_str = figfig.get_subplot((n_rows - 1, 0)).xaxis.get_label().get_text()

        # Set common y label
        y_lab_str = figfig.get_subplot((0, 0)).yaxis.get_label().get_text()
    else:
        x_lab_str = figfig.get_subplot().xaxis.get_label().get_text()

    # Replace x labels with one centered text
    for i in range(n_rows):
        for j in range(n_cols):
            figfig.get_subplot((i, j)).xaxis.get_label().set_text('')
    if grid_plot_pads is not None:
        figfig.fig.text(0.5, grid_plot_pads['bottom'], x_lab_str,
                        horizontalalignment='center',
                        verticalalignment='center')

    # Replace y labels with one centered text
    for i in range(n_rows):
        for j in range(n_cols):
            figfig.get_subplot((i, j)).yaxis.get_label().set_text('')
    if grid_plot_pads is not None:
        figfig.fig.text(grid_plot_pads['left'], 0.5, y_lab_str,
                        horizontalalignment='center',
                        verticalalignment='center',
                        rotation='vertical')

    # Add top title to each column to reflect varied parameter for each column
    for i in range(n_cols):
        ax = figfig.get_subplot((0, i))
        param_val = param_vals_x[i]
        if grid_plot_pads is not None:
            ax.text(0.5, 1 + grid_plot_pads['topMinor'], str(param_val),
                    horizontalalignment='center',
                    verticalalignment='bottom',
                    transform=ax.transAxes)
    if grid_plot_pads is not None and x_label_str_var is not None:
        figfig.fig.text(0.5, 1 - grid_plot_pads['top'], x_label_str_var,
                        horizontalalignment='center',
                        verticalalignment='center')

    # Add right title to each row to reflect varied parameter for each row
    if param_vals_y is not None:
        for i in range(n_rows):
            if param_vals_y[i] is None:
                continue
            ax = figfig.get_subplot((i, n_cols - 1))
            param_val = param_vals_y[i]
            if grid_plot_pads is not None:
                ax.text(1 + grid_plot_pads['rightMinor'], 0.5, str(param_val),
                        horizontalalignment='left',
                        verticalalignment='center',
                        transform=ax.transAxes,
                        rotation=270)
        if grid_plot_pads is not None and y_label_str_var is not None:
            figfig.fig.text(1 - grid_plot_pads['right'], 0.5, y_label_str_var,
                            horizontalalignment='center',
                            verticalalignment='center',
                            rotation=270)

    # Apply limits
    if x_lims is not None:
        [ax.set_xlim([x_lims[0], x_lims[1]]) for ax in figfig.axes.flat]
    if y_lims is not None:
        [ax.set_ylim([y_lims[0], y_lims[1]]) for ax in figfig.axes.flat]

    # Apply ticks
    if x_ticks is not None:
        [ax.set_xticks(x_ticks) for ax in figfig.axes.flat]
    if y_ticks is not None:
        [ax.set_yticks(y_ticks) for ax in figfig.axes.flat]

    # Set gridlines
    for ax in figfig.axes.flat:
        ax.grid(True, which='major', linewidth=1)
        ax.grid(True, which='minor', linewidth=0.1)
        ax.minorticks_on()
