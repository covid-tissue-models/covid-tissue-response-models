# todo - document stuff for easier usage by others (for now, use the workflow demo in CallableCoV2VTM.py!)

import os
import shutil
import csv
import matplotlib.pyplot as plt
import numpy as np

from PyQt5.QtCore import QObject

from cc3d.CompuCellSetup.CC3DCaller import CC3DCallerWorker
from cc3d.core import XMLUtils
from cc3d.core.BasicSimulationData import BasicSimulationData
from cc3d.core.GraphicsOffScreen.GenericDrawer import GenericDrawer
from cc3d.core.GraphicsUtils.ScreenshotManagerCore import ScreenshotManagerCore
from cc3d.cpp import PlayerPython
from cc3d.cpp.CompuCell import Dim3D
from cc3d.player5 import Configuration
from cc3d.player5.Simulation.CMLResultReader import CMLResultReader
from cc3d.player5.Utilities.utils import extract_address_int_from_vtk_object

export_data_desc = {'ir_data': ['ImmuneResp'],
                    'med_diff_data': ['MedViral',
                                      'MedCyt'],
                    'pop_data': ['Uninfected',
                                 'Infected',
                                 'InfectedSecreting',
                                 'Dying',
                                 'ImmuneCell',
                                 'ImmuneCellActivated'],
                    'spat_data': ['DeathComp',
                                  'InfectDist']}

x_label_str_transient = "Simulation time (MCS)"

y_label_str = {'ir_data': {'ImmuneResp': 'Immune response state variable'},
               'med_diff_data': {'MedViral': 'Total diffusive virus',
                                 'MedCyt': 'Total diffusive cytokine'},
               'pop_data': {'Uninfected': 'Number of uninfected cells',
                            'Infected': 'Number of infected cells',
                            'InfectedSecreting': 'Number of infected secreting cells',
                            'Dying': 'Number of dying cells',
                            'ImmuneCell': 'Number of immune cells',
                            'ImmuneCellActivated': 'Number of activated immune cells'},
               'spat_data': {'DeathComp': 'Cell death compactness (ul)',
                             'InfectDist': 'Infection distance (px)'}}

fig_save_names = {'ir_data': {'ImmuneResp': 'metric_immune_response_svar'},
                  'med_diff_data': {'MedViral': 'metric_diffusive_virus',
                                    'MedCyt': 'metric_diffusive_cytokine'},
                  'pop_data': {'Uninfected': 'metric_num_uninfected',
                               'Infected': 'metric_num_infected',
                               'InfectedSecreting': 'metric_num_infectedSecreting',
                               'Dying': 'metric_num_dying',
                               'ImmuneCell': 'metric_num_immune',
                               'ImmuneCellActivated': 'metric_num_immuneActivated'},
                  'spat_data': {'DeathComp': 'metric_death_compact',
                                'InfectDist': 'metric_infection_distance'}}


fig_suffix_trials = '_trials'
fig_suffix_stat = '_stat'


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
        sim_mcs = list(data_dict[list(data_dict.keys())[0]].keys())
        mean_data = dict()
        stdev_data = dict()
        for this_mcs in sim_mcs:
            mean_data[this_mcs] = dict()
            stdev_data[this_mcs] = dict()
            for param in param_names:
                this_data = []
                for trial_data in data_dict.values():
                    this_data.append(trial_data[this_mcs][param])

                mean_data[this_mcs][param] = float(np.average(this_data))
                stdev_data[this_mcs][param] = float(np.std(this_data))

        data_dict['batchMean'] = mean_data
        data_dict['batchStDev'] = stdev_data


def generate_batch_data_summary(cov2_vtm_sim_run):
    trial_dirs = [cov2_vtm_sim_run.get_run_output_dir(x) for x in range(cov2_vtm_sim_run.num_runs)]
    [convert_sim_data(trial_dir) for trial_dir in trial_dirs]
    batch_data_summary = {data_desc: collect_trial_data(data_desc, trial_dirs) for data_desc in export_data_desc.keys()}
    calculate_batch_data_stats(batch_data_summary)
    return batch_data_summary


def generate_transient_plot_trials(batch_data_summary, data_desc, var_name):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.grid()

    data_dict = batch_data_summary[data_desc]
    sim_mcs = list(data_dict[list(data_dict.keys())[0]].keys())
    for trial_idx in data_dict.keys():
        if trial_idx in ['batchMean', 'batchStDev']:
            continue
        y_data = [data_dict[trial_idx][this_mcs][var_name] for this_mcs in sim_mcs]
        ax.plot(sim_mcs, y_data, label='Trial {}'.format(trial_idx), marker='.')

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
    def __init__(self, cov2_vtm_sim_run):
        self.cov2_vtm_sim_run = cov2_vtm_sim_run

        # Structure is batch_data_summary[data_desc][trial_idx][this_mcs][param_name]
        self.batch_data_summary = generate_batch_data_summary(cov2_vtm_sim_run)
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

    def export_transient_plot_trials(self, loc=None):
        if loc is None:
            loc = self.cov2_vtm_sim_run.output_dir_root

        assert os.path.isdir(loc), "Results directory must be defined before rendering dump."

        fig_dir = self.get_fig_root_dir(loc, auto_make_dir=True)

        for data_desc in self.get_data_descs():
            for param_name in self.return_param_names(data_desc):
                fig, _ = self.generate_transient_plot_trials(data_desc, param_name)
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
    :param drawing_params:
    :return:
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


# todo - add parallel rendering
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

    def get_trial_vtk_dir(self, trial_idx):
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

            screenshot_spec_copy = os.path.join(ss_dir, os.path.basename(screenshot_spec))
            shutil.copyfile(screenshot_spec, screenshot_spec_copy)

    def load_trial_results(self, trial_idx):

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
        def get_screenshot_dir_name():
            return os.path.join(get_fig_spatial_dir(self.cov2_vtm_sim_run), f'run_{trial_idx}')

        self.scm.get_screenshot_dir_name = get_screenshot_dir_name
        self.scm.read_screenshot_description_file(ss_desc_file)

    def output_screenshots(self, mcs: int) -> None:
        self.scm.output_screenshots(mcs)

    def render_results(self):

        fig_spatial_dir = get_fig_spatial_dir(self.cov2_vtm_sim_run)
        fig_dir = os.path.dirname(fig_spatial_dir)

        if not os.path.isdir(fig_dir):
            os.mkdir(fig_dir)

        if not os.path.isdir(fig_spatial_dir):
            os.mkdir(fig_spatial_dir)

        trial_dirs = self.cov2_vtm_sim_run.get_trial_dirs()
        for trial_idx in range(len(trial_dirs)):
            print('CallableCC3DRenderer rendering trial {}'.format(trial_idx))
            self.load_trial_results(trial_idx)

            if self.cml_results_reader is None:
                print('No results loaded.')
                continue

            ss_dir = self.scm.get_screenshot_dir_name()
            if not os.path.isdir(ss_dir):
                os.mkdir(ss_dir)

            file_list = self.cml_results_reader.ldsFileList
            for file_number, file_name in enumerate(file_list):
                self.cml_results_reader.read_simulation_data_non_blocking(file_number)
                sim_data_int_addr = extract_address_int_from_vtk_object(self.cml_results_reader.simulationData)
                self.gd.field_extractor.setSimulationData(sim_data_int_addr)
                mcs = self.cml_results_reader.extract_mcs_number_from_file_name(file_name)
                print('...{}'.format(mcs))
                self.output_screenshots(mcs)
