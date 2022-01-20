# Written by T.J. Sego, Ph.D.

import sys
import os
from numpy import log10
from matplotlib.pyplot import FuncFormatter

sys.path.append(os.path.dirname(__file__))
import render_defs

# Export defs
rc_params = render_defs.rc_params_base.copy()
rc_params['figure.figsize'] = (3, 2.5)

export_fig_height = 2.5
export_fig_width = 3
params_export = ["Uninfected", "Infected", "InfectedSecreting", "Dying", "MedViral", "Viral"]
params_export = ['ImmuneResp', 'MedViral', 'Uninfected', "Infected", "InfectedSecreting", "Dying", 'r_max', 'Drug',
                 'ActiveMetabolite', 'Total_viral_RNA_in_cells', 'MedCyt', 'vir_auc', 'Mean_viral_RNA_in_cells']
# params_export = ['ImmuneResp', 'MedViral', 'Uninfected', "Infected", "InfectedSecreting", "Dying", 'r_max',
#                  'Total_viral_RNA_in_cells', 'MedCyt', 'vir_auc', 'Mean_viral_RNA_in_cells']
# params_export = ['ImmuneResp','Uninfected', "Infected", "InfectedSecreting", "Dying"]
#
# params_export = ['MedViral', 'Uninfected', "Infected", "InfectedSecreting", "Dying", 'r_max', 'Drug',
#                  'ActiveMetabolite', 'Total_viral_RNA_in_cells', 'MedCyt', 'vir_AUC', 'Mean_viral_RNA_in_cells']
# params_export = ['Drug',
#                  'ActiveMetabolite', 'Total_viral_RNA_in_cells', 'MedCyt', 'vir_AUC', 'Mean_viral_RNA_in_cells']

# params_export = ['Uninfected', "Infected", "InfectedSecreting", "Dying", 'r_max', 'Prodrug',
#                  'Metabolite4', 'tot_RNA']

s_to_mcs = 300.
exp_replicating_rate = 1.0 / 200.0 * 1.0 / 60.0
replicating_rate = exp_replicating_rate * s_to_mcs
# Plot manipulators
days_plot = [0, 14, 28]


# days_plot = [0, 28]


def hide_middle_tick(ax):
    # return
    for i, tick in enumerate(ax.get_xticklabels()):
        if i == 1:
            tick.set_visible(False)


def logs_tick_format_hide_middle(ax):
    ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    for i, tick in enumerate(ax.get_xticklabels()):
        if i == 1:
            tick.set_visible(False)


def exponent_ticks(value, tick_number):
    exp = log10(value)
    return int(exp)


def manip_plot_axes(fig, ax):
    render_defs.manip_plot_axis_conv(fig, ax)
    render_defs.manip_plot_time_label(fig, ax)

    xticks = [int(x) * 24 * 60 / 5 for x in days_plot]
    ax.set_xticks(xticks)
    # xticks = ax.xaxis.get_major_ticks()
    # xticks[1].label1.set_visible(False)
    # for i, tick in enumerate(ax.get_xticklabels()):
    #     if i == 1:
    #         tick.set_visible(False)


def manip_all(fig, ax):
    manip_plot_axes(fig, ax)
    render_defs.manip_plot_log(fig, ax)
    # render_defs.manip_plot_ode(fig, ax)


ec_plot_lim_manip = render_defs.plot_lim_manip(bottom=0, top=900)
ec_log_plot_lim_manip = render_defs.plot_lim_manip(bottom=1, top=500)
infected_lim_manip = render_defs.plot_lim_manip(bottom=0, top=500)
# ec_yticks = {'yticks': [1, 1E3/2, 1E3]}
ec_yticks = {'yticks': [0, 900]}

epit_lim_manip = render_defs.plot_lim_manip(bottom=0, top=1000)


def manip_susceptible(fig, ax):
    # epit_lim_manip(fig, ax)

    # manip_all(fig, ax)
    manip_plot_axes(fig, ax)
    ec_plot_lim_manip(fig, ax)
    hide_middle_tick(ax)
    render_defs.manip_set_ticks(**ec_yticks)(fig, ax)


def manip_infected(fig, ax):
    # manip_all(fig, ax)
    manip_plot_axes(fig, ax)
    infected_lim_manip(fig, ax)
    hide_middle_tick(ax)
    render_defs.manip_set_ticks(**ec_yticks)(fig, ax)


def manip_virus_releasing(fig, ax):
    # manip_all(fig, ax)
    manip_plot_axes(fig, ax)
    infected_lim_manip(fig, ax)
    hide_middle_tick(ax)
    render_defs.manip_set_ticks(**ec_yticks)(fig, ax)


def manip_dead(fig, ax):
    # manip_all(fig, ax)
    manip_plot_axes(fig, ax)
    ec_plot_lim_manip(fig, ax)
    hide_middle_tick(ax)
    render_defs.manip_set_ticks(**ec_yticks)(fig, ax)


def manip_immune_local(fig, ax):
    manip_plot_axes(fig, ax)
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=2)(fig, ax)


def manip_immune_lymph(fig, ax):
    manip_plot_axes(fig, ax)
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=2)(fig, ax)


def manip_virus(fig, ax):
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=1E6)(fig, ax)
    render_defs.manip_set_ticks(yticks=[1, 1E3, 1E6])(fig, ax)
    # ax.Axes.ticklabel_format(axis='y', style='sci')
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    logs_tick_format_hide_middle(ax)
    ax.tick_params(axis='y', labelsize=8)


def manip_cytokine(fig, ax):
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=1E8)(fig, ax)
    render_defs.manip_set_ticks(yticks=[1, 1E3, 1E7])(fig, ax)
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    logs_tick_format_hide_middle(ax)
    ax.tick_params(axis='y', labelsize=8)


def manip_cytokine_lymph(fig, ax):
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=2)(fig, ax)


def manip_death_virus(fig, ax):
    manip_all(fig, ax)
    hide_middle_tick(ax)
    ec_plot_lim_manip(fig, ax)


def manip_death_contact(fig, ax):
    manip_all(fig, ax)
    hide_middle_tick(ax)
    render_defs.plot_lim_manip(bottom=1, top=2)(fig, ax)


def manip_internal_RNA(fig, ax):
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=3000)(fig, ax)
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    logs_tick_format_hide_middle(ax)
    ax.tick_params(axis='y', labelsize=8)


def manip_mean_internal_RNA(fig, ax):
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=4)(fig, ax)
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    logs_tick_format_hide_middle(ax)
    ax.tick_params(axis='y', labelsize=8)


def manip_drug(fig, ax):
    manip_plot_axes(fig, ax)
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    logs_tick_format_hide_middle(ax)
    ax.tick_params(axis='y', labelsize=8)


def manip_rmax(fig, ax):
    # manip_plot_axes(fig, ax)
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1e-5, top=1.1 * replicating_rate)(fig, ax)
    hide_middle_tick(ax)
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    # render_defs.plot_lim_manip(bottom=1, top=0.025)(fig, ax)


def manip_vir_AUC(fig, ax):
    manip_all(fig, ax)
    render_defs.plot_lim_manip(bottom=1, top=35000)(fig, ax)
    # ax.yaxis.set_major_formatter(FuncFormatter(exponent_ticks))
    logs_tick_format_hide_middle(ax)
    ax.tick_params(axis='y', labelsize=8)

def manip_immune_svar(fig,ax):
    manip_plot_axes(fig, ax)
    render_defs.plot_lim_manip(bottom=-5, top=20)(fig, ax)
    hide_middle_tick(ax)



from collections import defaultdict

stat_plot_manips = defaultdict(lambda: None)
stat_plot_manips['Uninfected'] = manip_susceptible
stat_plot_manips['Infected'] = manip_infected
stat_plot_manips['InfectedSecreting'] = manip_virus_releasing
stat_plot_manips['Dying'] = manip_dead
stat_plot_manips['Immunelocal'] = manip_immune_local
stat_plot_manips['Immunelymphnode'] = manip_immune_lymph
stat_plot_manips['MedViral'] = manip_virus
stat_plot_manips['MedCyt'] = manip_cytokine
# stat_plot_manips['MedCytL'] = manip_cytokine_lymph
stat_plot_manips['Viral'] = manip_death_virus
stat_plot_manips['Contact'] = manip_death_contact
stat_plot_manips['Total_viral_RNA_in_cells'] = manip_internal_RNA
stat_plot_manips['Mean_viral_RNA_in_cells'] = manip_mean_internal_RNA
stat_plot_manips['r_max'] = manip_rmax
stat_plot_manips['Prodrug'] = manip_drug
stat_plot_manips['ActiveMetabolite'] = manip_drug
stat_plot_manips['vir_auc'] = manip_vir_AUC
stat_plot_manips['ImmuneResp'] = manip_immune_svar
