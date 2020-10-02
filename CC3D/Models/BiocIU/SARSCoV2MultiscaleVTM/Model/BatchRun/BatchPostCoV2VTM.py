# todo - document stuff for easier usage by others
import os
import sys
sys.path.append(os.environ['PYTHONPATH'])

import math
import shutil
import csv
try:
    import matplotlib.pyplot as plt
except ModuleNotFoundError:
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "matplotlib"])
    import matplotlib.pyplot as plt

import numpy as np
import vtk

from PyQt5.QtCore import QObject

from cc3d.core import XMLUtils
from cc3d.core.BasicSimulationData import BasicSimulationData
from cc3d.core.GraphicsOffScreen.GenericDrawer import GenericDrawer
from cc3d.core.GraphicsUtils.ScreenshotManagerCore import ScreenshotManagerCore
from cc3d.cpp import PlayerPython
from cc3d.cpp.CompuCell import Dim3D
from cc3d.player5 import Configuration
from cc3d.player5.Simulation.CMLResultReader import CMLResultReader
from cc3d.player5.Utilities.utils import extract_address_int_from_vtk_object
from Simulation.ViralInfectionVTMModelInputs import s_to_mcs

export_data_desc = {'ir_data': ['ImmuneResp'],
                    'med_diff_data': ['MedViral', 'MedCyt', 'MedOxi'],
                    'pop_data': ['Uninfected', 'Infected', 'InfectedSecreting', 'Dying', 'ImmuneCell',
                                 'ImmuneCellActivated'],
                    'spat_data': ['DeathComp', 'InfectDist'],
                    'death_data': ['Viral', 'OxiField', 'Contact', 'Bystander'],
                    'ddm_rmax_data': ['r_max'],
                    'ddm_data': ['Prodrug', 'Metabolite1', 'Metabolite2', 'Metabolite3', 'Metabolite4'],
                    'ddm_tot_RNA_data': ['tot_RNA']}

x_label_str_transient = "Simulation time (Days)"

y_label_str = {'ir_data': {'ImmuneResp': 'Immune response state variable'},
               'med_diff_data': {'MedViral': 'Total diffusive virus',
                                 'MedCyt': 'Total diffusive cytokine',
                                 'MedOxi': 'Total oxidative agent'},
               'pop_data': {'Uninfected': 'Number of uninfected cells',
                            'Infected': 'Number of infected cells',
                            'InfectedSecreting': 'Number of infected secreting cells',
                            'Dying': 'Number of dying cells',
                            'ImmuneCell': 'Number of immune cells',
                            'ImmuneCellActivated': 'Number of activated immune cells'},
               'spat_data': {'DeathComp': 'Cell death compactness (ul)',
                             'InfectDist': 'Infection distance (px)'},
               'death_data': {'Viral': 'Number of virally-induced apoptosis deaths',
                              'OxiField': 'Number of oxidative deaths',
                              'Contact': 'Number of cytotoxic kill deaths',
                              'Bystander': 'Number of bystander effect deaths'},
               'ddm_rmax_data': {'r_max': 'r_max value'},
               'ddm_data': {'Prodrug': 'Concentration of administered prodrug (A.U.)',
                            'Metabolite1': 'Concentration of 1st resulting metabolite (A.U.)',
                            'Metabolite2': 'Concentration of 2nd resulting metabolite (A.U.)',
                            'Metabolite3': 'Concentration of 3rd resulting metabolite (A.U.)',
                            'Metabolite4': 'Concentration of active metabolite (4th metabolite -- A.U.)'},
               'ddm_tot_RNA_data': {'tot_RNA': 'Total viral RNA in tissue cells'}
               }

fig_save_names = {'ir_data': {'ImmuneResp': 'metric_immune_response_svar'},
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
                  'ddm_data': {'Prodrug': 'metric_prodrug',
                               'Metabolite1': 'metric_metabolite1',
                               'Metabolite2': 'metric_metabolite2',
                               'Metabolite3': 'metric_metabolite3',
                               'Metabolite4': 'metric_metabolite4'},
                  'ddm_tot_RNA_data': {'tot_RNA': 'metric_total_RNA'}
                  }


fig_suffix_trials = '_trials'
fig_suffix_stat = '_stat'


def modify_data_desc(_file_root_name: str, _var_str_list: list) -> None:
    """
    Add or modify simulation data description in registry
    :param _file_root_name: root name of .dat file exported by simulation
    :param _var_str_list: list of strings of unique keys for each datum exported by simulation
    :return: None
    """
    export_data_desc[_file_root_name] = _var_str_list


def modify_y_label_str(_file_root_name: str, _var: str, _lab: str) -> None:
    """
    Add or modify y-axis labels in batch post-processing
    :param _file_root_name: root name of .dat file exported by simulation
    :param _var: key of datum
    :param _lab: new label
    :return: None
    """
    if _file_root_name not in y_label_str.keys():
        y_label_str[_file_root_name] = dict()
    y_label_str[_file_root_name][_var] = _lab


def modify_fig_save_names(_file_root_name: str, _var: str, _save_name) -> None:
    """
    Add or modify base name of saved figure for a datum in batch post-processing
    :param _file_root_name: root name of .dat file exported by simulation
    :param _var: key of datum
    :param _save_name: new figure base name
    :return: None
    """
    if _file_root_name not in fig_save_names:
        fig_save_names[_file_root_name] = dict()
    fig_save_names[_file_root_name][_var] = _save_name


def find_named_files(name: str, loc=os.path.dirname(__file__)) -> list:
    res_files = list()
    for root, dirs, names in os.walk(loc):
        [res_files.append(os.path.join(root, x)) for x in names if x == name]

    return res_files


def convert_files_2_csv(file_names: list) -> None:
    for name in file_names:
        name_csv = os.path.splitext(name)[0] + ".csv"
        os.rename(name, name_csv)


def convert_files_2_dat(file_names: list) -> None:
    for name in file_names:
        name_csv = os.path.splitext(name)[0] + ".dat"
        os.rename(name, name_csv)


def convert_sim_data(_loc):
    for _name in [x + '.dat' for x in export_data_desc.keys()]:
        convert_files_2_csv(find_named_files(_name, _loc))


def collect_trial_data(_export_name, _trial_dirs):
    trial_data = dict()
    param_names = export_data_desc[_export_name]
    num_trials = len(_trial_dirs)
    for trial_idx in range(num_trials):
        trial_data[trial_idx] = dict()
        trial_file = os.path.join(_trial_dirs[trial_idx], _export_name + '.csv')
        if not os.path.isfile(trial_file):
            trial_data[trial_idx] = None
            continue
        with open(trial_file) as csvfile:
            csv_reader = csv.reader(csvfile, delimiter=',')
            for row_data in csv_reader:
                this_mcs = int(row_data.pop(0))
                trial_data[trial_idx][this_mcs] = dict()
                for col_idx in range(len(row_data)):
                    trial_data[trial_idx][this_mcs][param_names[col_idx]] = float(row_data[col_idx])

    return trial_data


def calculate_batch_data_stats(batch_data_summary):
    for data_desc, data_dict in batch_data_summary.items():
        param_names = export_data_desc[data_desc]
        if data_dict[list(data_dict.keys())[0]] is None:
            continue
        sim_mcs = list(data_dict[list(data_dict.keys())[0]].keys())
        mean_data = dict()
        stdev_data = dict()
        for this_mcs in sim_mcs:
            mean_data[this_mcs] = dict()
            stdev_data[this_mcs] = dict()
            for param in param_names:
                this_data = []
                for trial_data in data_dict.values():
                    if trial_data is not None:
                        if this_mcs in trial_data.keys():
                            this_data.append(trial_data[this_mcs][param])

                if this_data:
                    mean_data[this_mcs][param] = float(np.average(this_data))
                    stdev_data[this_mcs][param] = float(np.std(this_data))

        data_dict['batchMean'] = mean_data
        data_dict['batchStDev'] = stdev_data


def generate_batch_data_summary(cov2_vtm_sim_run, step_list=None):
    trial_dirs = [cov2_vtm_sim_run.get_run_output_dir(x) for x in range(cov2_vtm_sim_run.num_runs)]
    [convert_sim_data(trial_dir) for trial_dir in trial_dirs]
    batch_data_summary = {data_desc: collect_trial_data(data_desc, trial_dirs) for data_desc in export_data_desc.keys()}
    # Filter data that wasn't found at all
    data_desc_found = {k: False for k in batch_data_summary.keys()}
    for data_desc, data_dict in batch_data_summary.items():
        for trial_dict in data_dict.values():
            if trial_dict is not None:
                data_desc_found[data_desc] = True
                continue
    batch_data_summary = {k: v for k, v in batch_data_summary.items() if data_desc_found[k]}
    # Apply step filter
    if step_list is not None:
        for data_desc, data_dict in batch_data_summary.items():
            for trial_idx, trial_dict in data_dict.items():
                if trial_dict is not None:
                    batch_data_summary[data_desc][trial_idx] = {k: v for k, v in trial_dict.items() if k in step_list}
    calculate_batch_data_stats(batch_data_summary)
    return batch_data_summary


def generate_transient_plot_trials(batch_data_summary, data_desc, var_name):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()

    data_dict = batch_data_summary[data_desc]
    # sim_mcs = list(data_dict[list(data_dict.keys())[0]].keys())
    for trial_idx in data_dict.keys():
        if trial_idx in ['batchMean', 'batchStDev']:
            continue
        sim_mcs = list(data_dict[trial_idx].keys())
        # sim_hours = float(sim_mcs[:]) * s_to_mcs / 60 / 60
        sim_days = [this_mcs*s_to_mcs/60/60/24 for this_mcs in sim_mcs]
        y_data = [data_dict[trial_idx][this_mcs][var_name] for this_mcs in sim_mcs]
        ax.plot(sim_days, y_data, label='Trial {}'.format(trial_idx), marker='.')
    ax.set_xlabel(x_label_str_transient)
    ax.set_ylabel(y_label_str[data_desc][var_name])
    # ax.legend()
    fig.tight_layout()

    return fig, ax


def generate_transient_plot_stat(batch_data_summary, data_desc, var_name, plot_stdev=True):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()

    data_dict = batch_data_summary[data_desc]
    sim_mcs = list(data_dict[list(data_dict.keys())[0]].keys())
    mean_dict = data_dict['batchMean']
    stdev_dict = data_dict['batchStDev']

    y_data = [mean_dict[this_mcs][var_name] for this_mcs in sim_mcs]
    ax.plot(sim_mcs, y_data, marker='.')
    if plot_stdev:
        yerr_data = [stdev_dict[this_mcs][var_name] for this_mcs in sim_mcs]
        ym_data = [y_data[i] - yerr_data[i] for i in range(len(y_data))]
        yp_data = [y_data[i] + yerr_data[i] for i in range(len(y_data))]
        ax.fill_between(sim_mcs, ym_data, yp_data, alpha=0.5)

    ax.set_xlabel(x_label_str_transient)
    ax.set_ylabel(y_label_str[data_desc][var_name])
    fig.tight_layout()

    return fig, ax


def generate_2var_plot_trials(batch_data_summary, var_name_hor, var_name_ver):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()

    data_desc_hor = None
    data_desc_ver = None
    for k, v in export_data_desc.items():
        if var_name_hor in v:
            data_desc_hor = k
        if var_name_ver in v:
            data_desc_ver = k

    assert data_desc_hor is not None, '{} is not a recognzied data variable'.format(data_desc_hor)
    assert data_desc_ver is not None, '{} is not a recognzied data variable'.format(data_desc_ver)

    data_dict_hor = batch_data_summary[data_desc_hor]
    data_dict_ver = batch_data_summary[data_desc_ver]
    sim_mcs = list(data_dict_hor[list(data_dict_hor.keys())[0]].keys())

    for trial_idx in data_dict_hor.keys():
        if trial_idx in ['batchMean', 'batchStDev']:
            continue
        x_data = [data_dict_hor[trial_idx][this_mcs][var_name_hor] for this_mcs in sim_mcs]
        y_data = [data_dict_ver[trial_idx][this_mcs][var_name_ver] for this_mcs in sim_mcs]
        ax.plot(x_data, y_data, label='Trial {}'.format(trial_idx), marker='.')

    ax.set_xlabel(y_label_str[data_desc_hor][var_name_hor])
    ax.set_ylabel(y_label_str[data_desc_ver][var_name_ver])
    # ax.legend()
    fig.tight_layout()

    return fig, ax


def generate_2var_plot_stat(batch_data_summary, var_name_hor, var_name_ver, plot_stdev=True):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()

    data_desc_hor = None
    data_desc_ver = None
    for k, v in export_data_desc.items():
        if var_name_hor in v:
            data_desc_hor = k
        if var_name_ver in v:
            data_desc_ver = k

    assert data_desc_hor is not None, '{} is not a recognzied data variable'.format(data_desc_hor)
    assert data_desc_ver is not None, '{} is not a recognzied data variable'.format(data_desc_ver)

    data_dict_hor = batch_data_summary[data_desc_hor]
    data_dict_ver = batch_data_summary[data_desc_ver]

    x_data = [data_dict_hor['batchMean'][k][var_name_hor]
              for k in batch_data_summary[data_desc_hor]['batchMean'].keys()]
    y_data = [data_dict_ver['batchMean'][k][var_name_ver]
              for k in batch_data_summary[data_desc_ver]['batchMean'].keys()]

    ax.plot(x_data, y_data, marker='.')
    if plot_stdev:
        xerr_data = [data_dict_hor['batchMean'][k][var_name_hor]
                     for k in batch_data_summary[data_desc_hor]['batchStDev'].keys()]
        yerr_data = [data_dict_ver['batchMean'][k][var_name_ver]
                     for k in batch_data_summary[data_desc_ver]['batchStDev'].keys()]

        xm_data = [x_data[i] - xerr_data[i] for i in range(len(xerr_data))]
        xp_data = [x_data[i] + xerr_data[i] for i in range(len(xerr_data))]
        ym_data = [y_data[i] - yerr_data[i] for i in range(len(yerr_data))]
        yp_data = [y_data[i] + yerr_data[i] for i in range(len(yerr_data))]
        ax.fill_between(x_data, ym_data, yp_data, alpha=0.5)
        ax.fill_betweenx(y_data, xm_data, xp_data, alpha=0.5, color='green')

    ax.set_xlabel(y_label_str[data_desc_hor][var_name_hor])
    ax.set_ylabel(y_label_str[data_desc_ver][var_name_ver])
    fig.tight_layout()

    return fig, ax


class CoV2VTMSimRunPost:
    """
    Renders simulation metrics data generated from executing a CallableCoV2VTM simulation batch
    """
    def __init__(self, cov2_vtm_sim_run, step_list=None):
        self.cov2_vtm_sim_run = cov2_vtm_sim_run

        # Structure is batch_data_summary[data_desc][trial_idx][this_mcs][param_name]
        self.batch_data_summary = generate_batch_data_summary(cov2_vtm_sim_run, step_list)
        self.__mcs_i = list(self.batch_data_summary[list(self.batch_data_summary.keys())[0]][0].keys())[0]

        self.__fig_suffix = '.png'

    def set_fig_suffix(self, _fig_suffix):
        self.__fig_suffix = _fig_suffix

    def get_data_descs(self):
        return list(self.batch_data_summary.keys())

    def return_param_names(self, data_desc):
        return list(self.batch_data_summary[data_desc][0][self.__mcs_i].keys())

    def get_fig_root_dir(self, _loc=None, auto_make_dir=False):
        if _loc is None:
            _loc = self.cov2_vtm_sim_run.output_dir_root

        fig_dir = os.path.join(_loc, f'Figs', f'Metrics')

        if auto_make_dir:
            if not os.path.isdir(fig_dir):
                if not os.path.isdir(os.path.join(_loc, f'Figs')):
                    os.mkdir(os.path.join(_loc, f'Figs'))
                os.mkdir(fig_dir)

        assert os.path.isdir(fig_dir)

        return fig_dir

    def generate_transient_plot_trials(self, data_desc, param_name):
        fig, ax = generate_transient_plot_trials(self.batch_data_summary, data_desc, param_name)
        return fig, ax

    def generate_transient_plot_stat(self, data_desc, param_name, plot_stdev=True):
        fig, ax = generate_transient_plot_stat(self.batch_data_summary, data_desc, param_name, plot_stdev)
        return fig, ax

    def generate_2var_plot_trials(self, var_name_hor, var_name_ver):
        fig, ax = generate_2var_plot_trials(self.batch_data_summary, var_name_hor, var_name_ver)
        return fig, ax

    def generate_2var_plot_stat(self, var_name_hor, var_name_ver, plot_stdev=True):
        fig, ax = generate_2var_plot_stat(self.batch_data_summary, var_name_hor, var_name_ver, plot_stdev)
        return fig, ax

    def generate_transient_plot_trials_filename(self, data_desc, param_name, fig_suffix=None, fig_dir=None):
        if fig_dir is None:
            fig_dir = self.get_fig_root_dir()

        if fig_suffix is None:
            fig_suffix = self.__fig_suffix

        fig_save_name_rel = fig_save_names[data_desc][param_name] + fig_suffix_trials + fig_suffix
        return os.path.join(fig_dir, fig_save_name_rel)

    def generate_transient_plot_stat_filename(self, data_desc, param_name, fig_suffix=None, fig_dir=None):
        if fig_dir is None:
            fig_dir = self.get_fig_root_dir()

        if fig_suffix is None:
            fig_suffix = self.__fig_suffix

        fig_save_name_rel = fig_save_names[data_desc][param_name] + fig_suffix_stat + fig_suffix
        return os.path.join(fig_dir, fig_save_name_rel)

    def generate_2var_plot_trials_filename(self, var_name_hor, var_name_ver, fig_suffix=None, fig_dir=None):
        if fig_dir is None:
            fig_dir = self.get_fig_root_dir()

        if fig_suffix is None:
            fig_suffix = self.__fig_suffix

        fig_save_name_rel = 'metric_' + var_name_hor + '_and_' + var_name_ver + fig_suffix_trials + fig_suffix
        return os.path.join(fig_dir, fig_save_name_rel)

    def generate_2var_plot_stat_filename(self, var_name_hor, var_name_ver, fig_suffix=None, fig_dir=None):
        if fig_dir is None:
            fig_dir = self.get_fig_root_dir()

        if fig_suffix is None:
            fig_suffix = self.__fig_suffix

        fig_save_name_rel = 'metric_' + var_name_hor + '_and_' + var_name_ver + fig_suffix_stat + fig_suffix
        return os.path.join(fig_dir, fig_save_name_rel)

    def export_transient_plot_trials(self, loc=None, manipulators=None):
        if loc is None:
            loc = self.cov2_vtm_sim_run.output_dir_root

        assert os.path.isdir(loc), "Results directory must be defined before rendering dump."
        assert manipulators is None or isinstance(manipulators,
                                                  dict), "manipulators must be None or a dictionary of manipulator functions"

        fig_dir = self.get_fig_root_dir(loc, auto_make_dir=True)

        for data_desc in self.get_data_descs():
            for param_name in self.return_param_names(data_desc):
                fig, ax = self.generate_transient_plot_trials(data_desc, param_name)
                if manipulators is not None and param_name in manipulators.keys():
                    manipulators[param_name](fig, ax)
                fig.savefig(self.generate_transient_plot_trials_filename(data_desc, param_name, fig_dir=fig_dir))
                plt.close(fig)

    def export_transient_plot_stat(self, loc=None, plot_stdev=True):
        if loc is None:
            loc = self.cov2_vtm_sim_run.output_dir_root

        assert os.path.isdir(loc), "Results directory must be defined before rendering dump."

        fig_dir = self.get_fig_root_dir(loc, auto_make_dir=True)

        for data_desc in self.get_data_descs():
            for param_name in self.return_param_names(data_desc):
                fig, _ = self.generate_transient_plot_stat(data_desc, param_name, plot_stdev)
                fig.savefig(self.generate_transient_plot_stat_filename(data_desc, param_name, fig_dir=fig_dir))
                plt.close(fig)

    def export_2var_plot_trials(self, var_name_hor, var_name_ver, loc=None):
        if loc is None:
            loc = self.cov2_vtm_sim_run.output_dir_root

        assert os.path.isdir(loc), "Results directory must be defined before rendering dump."

        fig_dir = self.get_fig_root_dir(loc, auto_make_dir=True)

        fig, _ = self.generate_2var_plot_trials(var_name_hor, var_name_ver)
        fig.savefig(self.generate_2var_plot_trials_filename(var_name_hor, var_name_ver, fig_dir=fig_dir))
        plt.close(fig)

    def export_2var_plot_stat(self, var_name_hor, var_name_ver, loc=None, plot_stdev=True):
        if loc is None:
            loc = self.cov2_vtm_sim_run.output_dir_root

        assert os.path.isdir(loc), "Results directory must be defined before rendering dump."

        fig_dir = self.get_fig_root_dir(loc, auto_make_dir=True)

        fig, _ = self.generate_2var_plot_stat(var_name_hor, var_name_ver, plot_stdev)
        fig.savefig(self.generate_2var_plot_stat_filename(var_name_hor, var_name_ver, fig_dir=fig_dir))
        plt.close(fig)


class CC3DUIDummy(QObject):
    """
    Some trickery to fake the launching of Player
    """
    def __init__(self, field_dim: Dim3D):
        super().__init__()
        self.fieldExtractor = PlayerPython.FieldExtractorCML()
        self.fieldExtractor.setFieldDim(field_dim)


def get_lattice_description_file(loc):
    """
    Generate CML Results Reader instance for a results directory
    :param loc: Path to directory with lattice description file and results
    :return: path to lattice description file
    """
    lds_file = None
    for file in os.listdir(loc):
        if file.endswith('.dml'):
            lds_file = os.path.join(loc, file)
            break

    if lds_file is None:
        'Could not locate lattice description file in ' + loc

    return lds_file


def get_results_reader_no_ui(lds_file):
    """
    Generate CML Results Reader instance for a results lattice description file
    :param lds_file: Path to lattice description file
    :return: CMLResultsReader instance
    """
    xml2_obj_converter = XMLUtils.Xml2Obj()
    root_element = xml2_obj_converter.Parse(lds_file)
    dim_element = root_element.getFirstElement("Dimensions")
    field_dim = Dim3D()
    field_dim.x = int(dim_element.getAttribute("x"))
    field_dim.y = int(dim_element.getAttribute("y"))
    field_dim.z = int(dim_element.getAttribute("z"))

    ui_dummy = CC3DUIDummy(field_dim)

    return CMLResultReader(ui_dummy), ui_dummy


def get_fig_spatial_dir(cov2_vtm_sim_run):
    return os.path.join(cov2_vtm_sim_run.output_dir_root, f'Figs', f'Spatial')


def get_trial_vtk_dir(cov2_vtm_sim_run, trial_idx):
    return os.path.join(cov2_vtm_sim_run.get_run_output_dir(trial_idx), f'LatticeData')


def get_trial_vtk_mcs_list(cov2_vtm_sim_run, trial_idx):
    vtk_files = [x for x in os.listdir(get_trial_vtk_dir(cov2_vtm_sim_run, trial_idx)) if x.endswith('.vtk')]
    return [int(x.rstrip('.vtk').split('_')[-1]) for x in vtk_files]


class GenericDrawerFree(GenericDrawer):
    """
    Removes dependency on persistent globals
    """
    def __init__(self, parent=None, originating_widget=None):
        super().__init__(parent, originating_widget)

    def get_model_view(self, drawing_params):
        model, view = GenericDrawer.get_model_view(self, drawing_params)

        lattice_type = self.lattice_type
        lattice_type_str = [k for k, v in Configuration.LATTICE_TYPES.items() if v == lattice_type][0]

        def init_lattice_type():
            model.lattice_type = lattice_type
            model.lattice_type_str = lattice_type_str

        model.init_lattice_type = init_lattice_type
        model.init_lattice_type()

        return model, view


class CallableCC3DRenderer:
    """
    Performs CC3D rendering of data generated from executing a CallableCoV2VTM simulation batch without launching Player
    """
    def __init__(self, cov2_vtm_sim_run):
        self.cov2_vtm_sim_run = cov2_vtm_sim_run

        self.gd = GenericDrawerFree()
        self.scm = ScreenshotManagerCore()
        self.scm.gd = self.gd

        self.cml_results_reader = None

        # Methods for modifying specification of GenericDrawer
        self._gd_manipulators = {}

        # Methods for modifying specification of ScreenshotData
        self._sc_manipulators = {}

    def get_trial_vtk_dir(self, trial_idx):
        """
        Returns path to directory where exported vtk files from simulation should be found
        :param trial_idx: index of a trial
        :return: path to directory containing exported vtk files from simulation
        """
        return get_trial_vtk_dir(self.cov2_vtm_sim_run, trial_idx)

    def load_screenshot_specs(self, screenshot_spec, trial_idx=None):
        """
        Loads screenshot specifications for rendering
        :param screenshot_spec: path to json screenshot specification (can be generated in Player)
        :param trial_idx: trial for which to apply the specification; default is all runs of the loaded batch
        :return: None
        """
        if trial_idx is None:
            trial_vtk_dirs = [self.get_trial_vtk_dir(i) for i in range(self.cov2_vtm_sim_run.num_runs)]
        else:
            trial_vtk_dirs = [self.get_trial_vtk_dir(trial_idx)]

        screenshot_spec = os.path.abspath(screenshot_spec)
        for trial_vtk_dir in trial_vtk_dirs:
            ss_dir = os.path.join(trial_vtk_dir, 'screenshot_data')
            if not os.path.isdir(ss_dir):
                os.mkdir(ss_dir)

            screenshot_spec_copy = os.path.join(ss_dir, 'screenshots.json')
            shutil.copyfile(screenshot_spec, screenshot_spec_copy)

    def load_trial_results(self, trial_idx):
        """
        Load results for a trial of a batch into memory; must be executed before manipulating rendering specs
        :param trial_idx: index of trial
        :return: None
        """

        lds_loc = self.get_trial_vtk_dir(trial_idx)

        lds_file = get_lattice_description_file(lds_loc)

        if lds_file is None:
            return

        self.cml_results_reader, ui_dummy = get_results_reader_no_ui(lds_file)

        if self.cml_results_reader is None:
            return

        self.gd.set_field_extractor(ui_dummy.fieldExtractor)
        self.cml_results_reader.extract_lattice_description_info(lds_file)

        ss_desc_file = None
        for root, dirs, names in os.walk(lds_loc):
            for name in names:
                if name == 'screenshots.json':
                    ss_desc_file = os.path.join(root, name)
                    break

        if ss_desc_file is None:
            print('No screenshot description found.')
            return

        self.scm.bsd = BasicSimulationData()
        self.scm.bsd.fieldDim = self.cml_results_reader.fieldDim
        self.scm.bsd.numberOfSteps = self.cml_results_reader.numberOfSteps

        # Overload static ScreenshotManagerCore.get_screenshot_dir_name, since it relies on persistent_globals
        screenshot_dir_name = os.path.join(get_fig_spatial_dir(self.cov2_vtm_sim_run), f'run_{trial_idx}')

        def get_screenshot_dir_name():
            return screenshot_dir_name

        self.scm.get_screenshot_dir_name = get_screenshot_dir_name
        self.scm.read_screenshot_description_file(ss_desc_file)

    def output_screenshots(self, mcs: int) -> None:
        """
        Executes screenshot rendering and write to disk for a simulation step of a trial in memory
        :param mcs: simulation step
        :return: None
        """
        self.scm.output_screenshots(mcs)

    def prep_output_dir(self):
        """
        Prep directory for output from rendering
        :return: None
        """
        fig_spatial_dir = get_fig_spatial_dir(self.cov2_vtm_sim_run)
        fig_dir = os.path.dirname(fig_spatial_dir)

        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

        if not os.path.isdir(fig_spatial_dir):
            os.mkdir(fig_spatial_dir)

    # todo - add API for defining rendering specs so users don't have to search through the details of GenericDrawer;
    #  API should include convenience function for retrieving current specs
    def load_rendering_manipulator(self, gd_manipulator, trial_idx=0, mcs=0):
        """
        Loads a function for manipulating GenericDrawer rendering parameters at a simulation step of a trial
        Manipulations are applied to all subsequent rendering processes
        :param gd_manipulator: manipulator; signature should have argument of GenericDrawer object
        :param trial_idx: trial at which to apply the manipulator
        :param mcs: step at which to apply the manipulator
        :return: None
        """
        if trial_idx not in self._gd_manipulators.keys():
            self._gd_manipulators[trial_idx] = dict()

        self._gd_manipulators[trial_idx][mcs] = gd_manipulator

    # todo - add API for defining rendering specs so users don't have to search through the details of GenericDrawer;
    #  API should include convenience function for retrieving current specs
    def load_screenshot_manipulator(self, sc_manipulator, trial_idx=0, mcs=0):
        """
        Loads a function for manipulating ScreenshotManagerCore rendering parameters at a simulation step of a trial
        Manipulations are applied to all subsequent rendering processes
        :param sc_manipulator: manipulator; signature should have argument of ScreenshotManagerCore object
        :param trial_idx: trial at which to apply the manipulator
        :param mcs: step at which to apply the manipulator
        :return: None
        """
        if trial_idx not in self._sc_manipulators.keys():
            self._sc_manipulators[trial_idx] = dict()

        self._sc_manipulators[trial_idx][mcs] = sc_manipulator

    def get_results_min_max(self, trial_idx):
        """
        Gets minimum and maximum over all simulation time for all available results
        :return: {dict} range per available field
        """
        min_max_dict = dict()

        self.load_trial_results(trial_idx)

        if self.cml_results_reader is None:
            print('No results loaded for trial {}.'.format(trial_idx))
            return None

        file_list = self.cml_results_reader.ldsFileList
        for file_number, file_name in enumerate(file_list):
            self.cml_results_reader.read_simulation_data_non_blocking(file_number)
            sim_data_int_addr = extract_address_int_from_vtk_object(self.cml_results_reader.simulationData)
            self.gd.field_extractor.setSimulationData(sim_data_int_addr)

            for field_name, screenshot_data in self.scm.screenshotDataDict.items():
                min_max = self._get_field_min_max(screenshot_data)
                if min_max is not None:
                    if field_name not in min_max_dict.keys():
                        min_max_dict[field_name] = min_max
                    else:
                        if min_max[0] < min_max_dict[field_name][0]:
                            min_max_dict[field_name][0] = min_max[0]
                        if min_max[1] > min_max_dict[field_name][1]:
                            min_max_dict[field_name][1] = min_max[1]

        return min_max_dict

    def _get_field_min_max(self, screenshot_data):
        """
        Gets minimum and maximum of field described in screenshot data; returns None if unavailable
        :param screenshot_data: screenshot data for a field
        :return: minimum and maximum
        """

        field_name = screenshot_data.plotData[0]
        fieldType = screenshot_data.plotData[1]
        plane = screenshot_data.projection
        planePos = screenshot_data.projectionPosition

        con_array = vtk.vtkDoubleArray()
        con_array.SetName("concentration")
        con_array_int_addr = extract_address_int_from_vtk_object(vtkObj=con_array)

        field_type = fieldType.lower()
        if field_type == 'confield':
            fill_successful = self.gd.field_extractor.fillConFieldData2D(con_array_int_addr,
                                                                         field_name,
                                                                         plane,
                                                                         planePos)
        elif field_type == 'scalarfield':
            fill_successful = self.gd.field_extractor.fillScalarFieldData2D(con_array_int_addr,
                                                                            field_name,
                                                                            plane,
                                                                            planePos)
        elif field_type == 'scalarfieldcelllevel':
            fill_successful = self.gd.field_extractor.fillScalarFieldCellLevelData2D(con_array_int_addr,
                                                                                     field_name,
                                                                                     plane,
                                                                                     planePos)

        else:
            return None

        if not fill_successful:
            return None

        con_array_range = con_array.GetRange()
        min_max = [con_array_range[0], con_array_range[1]]
        return min_max

    def _get_rendering_manipulator(self, trial_idx, mcs):
        if trial_idx not in self._gd_manipulators.keys() or mcs not in self._gd_manipulators[trial_idx].keys():
            return None
        else:
            return self._gd_manipulators[trial_idx][mcs]

    def _get_screenshot_manipulator(self, trial_idx, mcs):
        if trial_idx not in self._sc_manipulators.keys() or mcs not in self._sc_manipulators[trial_idx].keys():
            return None
        else:
            return self._sc_manipulators[trial_idx][mcs]

    def _render_trial(self, trial_idx):
        """
        Main routine to perform rendering for a trial from batch run
        :param trial_idx: index of trial
        :return: None
        """
        print('CallableCC3DRenderer rendering trial {}'.format(trial_idx))
        self.load_trial_results(trial_idx)

        if self.cml_results_reader is None:
            print('No results loaded.')
            return

        ss_dir = self.scm.get_screenshot_dir_name()
        if not os.path.isdir(ss_dir):
            os.mkdir(ss_dir)

        file_list = self.cml_results_reader.ldsFileList
        for file_number, file_name in enumerate(file_list):
            self.cml_results_reader.read_simulation_data_non_blocking(file_number)
            sim_data_int_addr = extract_address_int_from_vtk_object(self.cml_results_reader.simulationData)
            mcs = self.cml_results_reader.extract_mcs_number_from_file_name(file_name)
            self.gd.field_extractor.setSimulationData(sim_data_int_addr)
            gd_manipulator = self._get_rendering_manipulator(trial_idx, mcs)
            if gd_manipulator is not None:
                gd_manipulator(self.gd)
            print('...{}'.format(mcs))
            self.output_screenshots(mcs)

    def render_results(self):
        """
        Render all trials from batch run
        :return: None
        """
        self.prep_output_dir()
        [self._render_trial(trial_idx) for trial_idx in range(len(self.cov2_vtm_sim_run.get_trial_dirs()))]

    def render_trial_results(self, trial_idx):
        """
        Render a trial from batch run
        :param trial_idx: index of trial
        :return: None
        """
        self.prep_output_dir()
        self._render_trial(trial_idx)

    def render_trial_results_par(self, opts=None):
        """
        Render all trials in parallel
        :return: None
        """
        import multiprocessing
        from cc3d.CompuCellSetup.CC3DCaller import CC3DCallerWorker

        # Start workers
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        workers = [CC3DCallerWorker(tasks, results) for _ in range(self.cov2_vtm_sim_run.num_workers)]
        [w.start() for w in workers]

        # Enqueue jobs
        num_runs = len(self.cov2_vtm_sim_run.get_trial_dirs())
        [tasks.put(_RenderJob(self.cov2_vtm_sim_run, run_idx, opts)) for run_idx in range(num_runs)]

        # Add a stop task for each of worker
        [tasks.put(None) for _ in workers]

        tasks.join()


class _RenderJob:
    def __init__(self, _cov2_vtm_sim_run, _run_idx, opts=None):
        super().__init__()
        self._cov2_vtm_sim_run = _cov2_vtm_sim_run
        self._run_idx = _run_idx
        if opts is not None:
            self._opts = opts
        else:
            self._opts = dict()

    def run(self):
        print(f'Rendering job: Run {self._run_idx}')
        try:
            _renderer = CallableCC3DRenderer(self._cov2_vtm_sim_run)
            _renderer.load_trial_results(self._run_idx)

            if 'log_scale' in self._opts.keys() and self._opts['log_scale']:
                # Apply log scale to all field renders
                def gd_manipulator(gd):
                    gd.draw_model_2D.clut.SetScaleToLog10()
                _renderer.load_rendering_manipulator(gd_manipulator)

            if 'fixed_caxes' in self._opts.keys() and self._opts['fixed_caxes']:
                # Apply fixed legends to all field renders
                min_max_dict = _renderer.get_results_min_max(self._run_idx)

                def sc_manipulator(scm):
                    for field_name, min_max in min_max_dict.items():
                        scm.screenshotDataDict[field_name].metadata['MinRangeFixed'] = True
                        scm.screenshotDataDict[field_name].metadata['MinRange'] = math.ceil(
                            min_max[1] * 10) / 10 / 1E6
                        scm.screenshotDataDict[field_name].metadata['MaxRangeFixed'] = True
                        scm.screenshotDataDict[field_name].metadata['MaxRange'] = math.ceil(min_max[1] * 10) / 10

                _renderer.load_screenshot_manipulator(sc_manipulator)

            _renderer.render_trial_results(self._run_idx)
            return True
        except Exception:
            return False


class CallableCC3DDataRenderer(CallableCC3DRenderer):
    """
    Performs CC3D rendering of data generated from executing a CallableCoV2VTM simulation batch without launching Player
    Like CallableCC3DRenderer, but works on individual directories of data instead of a CallableCoV2VTM instance
    """
    def __init__(self, data_dirs, out_dirs, set_labs=None, run_labs=None, num_workers=1):
        super().__init__(None)

        self.data_dirs = data_dirs
        self.out_dirs = out_dirs
        if set_labs is not None:
            self.set_labs = set_labs
        else:
            self.set_labs = [0] * len(self.data_dirs)
        if run_labs is not None:
            self.run_labs = run_labs
        else:
            self.run_labs = [0] * len(self.data_dirs)
        self.num_workers = num_workers

    def get_trial_vtk_dir(self, trial_idx):
        """
        Returns path to directory where exported vtk files from simulation should be found
        :param trial_idx: index of a trial
        :return: path to directory containing exported vtk files from simulation
        """
        return os.path.join(self.data_dirs[trial_idx], 'LatticeData')

    def get_fig_spatial_dir(self, trial_idx):
        return os.path.join(self.out_dirs[trial_idx], f'set_{self.set_labs[trial_idx]}', 'Figs', 'Spatial')

    def load_screenshot_specs(self, screenshot_spec, trial_idx=None):
        """
        Loads screenshot specifications for rendering
        :param screenshot_spec: path to json screenshot specification (can be generated in Player)
        :param trial_idx: trial for which to apply the specification; default is all runs of the loaded batch
        :return: None
        """
        if trial_idx is None:
            trial_vtk_dirs = [self.get_trial_vtk_dir(i) for i in range(len(self.data_dirs))]
        else:
            trial_vtk_dirs = [self.get_trial_vtk_dir(trial_idx)]

        screenshot_spec = os.path.abspath(screenshot_spec)
        for trial_vtk_dir in trial_vtk_dirs:
            ss_dir = os.path.join(trial_vtk_dir, 'screenshot_data')
            if not os.path.isdir(ss_dir):
                os.mkdir(ss_dir)

            screenshot_spec_copy = os.path.join(ss_dir, 'screenshots.json')
            shutil.copyfile(screenshot_spec, screenshot_spec_copy)

    def load_trial_results(self, trial_idx):
        """
        Load results for a trial of a batch into memory; must be executed before manipulating rendering specs
        :param trial_idx: index of trial
        :return: None
        """

        lds_loc = self.get_trial_vtk_dir(trial_idx)

        lds_file = get_lattice_description_file(lds_loc)

        if lds_file is None:
            return

        self.cml_results_reader, ui_dummy = get_results_reader_no_ui(lds_file)

        if self.cml_results_reader is None:
            return

        self.gd.set_field_extractor(ui_dummy.fieldExtractor)
        self.cml_results_reader.extract_lattice_description_info(lds_file)

        ss_desc_file = None
        for root, dirs, names in os.walk(lds_loc):
            for name in names:
                if name == 'screenshots.json':
                    ss_desc_file = os.path.join(root, name)
                    break

        if ss_desc_file is None:
            print('No screenshot description found.')
            return

        self.scm.bsd = BasicSimulationData()
        self.scm.bsd.fieldDim = self.cml_results_reader.fieldDim
        self.scm.bsd.numberOfSteps = self.cml_results_reader.numberOfSteps

        # Overload static ScreenshotManagerCore.get_screenshot_dir_name, since it relies on persistent_globals
        screenshot_dir_name = os.path.join(self.get_fig_spatial_dir(trial_idx), f'run_{self.run_labs[trial_idx]}')

        def get_screenshot_dir_name():
            return screenshot_dir_name

        self.scm.get_screenshot_dir_name = get_screenshot_dir_name
        self.scm.read_screenshot_description_file(ss_desc_file)

    def prep_output_dir(self):
        """
        Prep directory for output from rendering
        :return: None
        """
        [os.makedirs(self.get_fig_spatial_dir(t))
         for t in range(len(self.out_dirs)) if not os.path.isdir(self.get_fig_spatial_dir(t))]

    def get_results_min_max(self, trial_idx):
        """
        Gets minimum and maximum over all simulation time for all available results
        :return: {dict} range per available field
        """
        min_max_dict = dict()

        self.load_trial_results(trial_idx)

        if self.cml_results_reader is None:
            print('No results loaded for trial {}.'.format(trial_idx))
            return None

        file_list = self.cml_results_reader.ldsFileList
        for file_number, file_name in enumerate(file_list):
            self.cml_results_reader.read_simulation_data_non_blocking(file_number)
            sim_data_int_addr = extract_address_int_from_vtk_object(self.cml_results_reader.simulationData)
            self.gd.field_extractor.setSimulationData(sim_data_int_addr)

            for field_name, screenshot_data in self.scm.screenshotDataDict.items():
                min_max = self._get_field_min_max(screenshot_data)
                if min_max is not None:
                    if field_name not in min_max_dict.keys():
                        min_max_dict[field_name] = min_max
                    else:
                        if min_max[0] < min_max_dict[field_name][0]:
                            min_max_dict[field_name][0] = min_max[0]
                        if min_max[1] > min_max_dict[field_name][1]:
                            min_max_dict[field_name][1] = min_max[1]

        return min_max_dict

    def _render_trial(self, trial_idx):
        """
        Main routine to perform rendering for a trial from batch run
        :param trial_idx: index of trial
        :return: None
        """
        print('CallableCC3DRenderer rendering trial {}'.format(trial_idx))
        self.load_trial_results(trial_idx)

        if self.cml_results_reader is None:
            print('No results loaded.')
            return

        ss_dir = self.scm.get_screenshot_dir_name()
        if not os.path.isdir(ss_dir):
            os.mkdir(ss_dir)

        file_list = self.cml_results_reader.ldsFileList
        for file_number, file_name in enumerate(file_list):
            self.cml_results_reader.read_simulation_data_non_blocking(file_number)
            sim_data_int_addr = extract_address_int_from_vtk_object(self.cml_results_reader.simulationData)
            mcs = self.cml_results_reader.extract_mcs_number_from_file_name(file_name)
            self.gd.field_extractor.setSimulationData(sim_data_int_addr)
            gd_manipulator = self._get_rendering_manipulator(trial_idx, mcs)
            if gd_manipulator is not None:
                gd_manipulator(self.gd)
            print('...{}'.format(mcs))
            self.output_screenshots(mcs)

    def render_results(self):
        """
        Render all trials from batch run
        :return: None
        """
        self.prep_output_dir()
        [self._render_trial(trial_idx) for trial_idx in range(len(self.data_dirs))]

    def render_trial_results(self, trial_idx):
        """
        Render a trial from batch run
        :param trial_idx: index of trial
        :return: None
        """
        self.prep_output_dir()
        self._render_trial(trial_idx)

    def render_trial_results_par(self, opts=None):
        """
        Render all trials in parallel
        :return: None
        """
        import multiprocessing
        from cc3d.CompuCellSetup.CC3DCaller import CC3DCallerWorker

        # Start workers
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        workers = [CC3DCallerWorker(tasks, results) for _ in range(self.num_workers)]
        [w.start() for w in workers]

        # Enqueue jobs
        import time
        for r in range(len(self.data_dirs)):
            while tasks.full():
                time.sleep(1)
            tasks.put(_RenderDataJob(self.data_dirs[r], self.out_dirs[r], self.set_labs[r], self.run_labs[r], opts))

        # Add a stop task for each of worker
        for _ in workers:
            while tasks.full():
                time.sleep(1)
            tasks.put(None)

        tasks.join()


class _RenderDataJob:
    def __init__(self, _data_dir, _out_dir, _set_lab, _run_lab, opts=None):
        self._data_dir = _data_dir
        self._out_dir = _out_dir
        self._set_lab = _set_lab
        self._run_lab = _run_lab
        if opts is not None:
            self._opts = opts
        else:
            self._opts = dict()

    def run(self):
        try:
            _renderer = CallableCC3DDataRenderer(data_dirs=[self._data_dir],
                                                 out_dirs=[self._out_dir],
                                                 set_labs=[self._set_lab],
                                                 run_labs=[self._run_lab])
            _renderer.load_trial_results(0)

            if 'log_scale' in self._opts.keys() and self._opts['log_scale']:
                # Apply log scale to all field renders
                def gd_manipulator(gd):
                    gd.draw_model_2D.clut.SetScaleToLog10()
                _renderer.load_rendering_manipulator(gd_manipulator)

            if 'fixed_caxes' in self._opts.keys() and self._opts['fixed_caxes']:
                # Apply fixed legends to all field renders
                min_max_dict = _renderer.get_results_min_max(0)

                def sc_manipulator(scm):
                    for field_name, min_max in min_max_dict.items():
                        scm.screenshotDataDict[field_name].metadata['MinRangeFixed'] = True
                        scm.screenshotDataDict[field_name].metadata['MinRange'] = math.ceil(
                            min_max[1] * 10) / 10 / 1E6
                        scm.screenshotDataDict[field_name].metadata['MaxRangeFixed'] = True
                        scm.screenshotDataDict[field_name].metadata['MaxRange'] = math.ceil(min_max[1] * 10) / 10

                _renderer.load_screenshot_manipulator(sc_manipulator)

            _renderer.render_trial_results(0)
            return True
        except Exception:
            return False
