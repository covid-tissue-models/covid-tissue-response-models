import multiprocessing
import shutil
import sys

# Key where Python API input dictionary is stored
cc3d_input_key = '__input_dict__'

# Key in Python API input dictionary where batch configuration inputs are stored for batch runs
cc3d_batch_key = '__cc3d_batch__'

# Key for sharing auto inputs over multiple processes
cc3d_auto_key = '__auto_inputs__'

# Automatic tasks

#   Plot variables by module
#       These are automatically disabled when running in batch mode, so the user doesn't have to modify any simulation
#       scripts when running a batch
mod_plot_vars = {'ViralInfectionVTMModelInputs': ['plot_vrm_data_freq', 'plot_vrm_data_freq', 'plot_vim_data_freq',
                                                  'plot_pop_data_freq', 'plot_ir_data_freq', 'plot_med_diff_data_freq',
                                                  'plot_spat_data_freq', 'plot_death_data_freq'],
                 'Models.DrugDosingModel.DrugDosingInputs': ['plot_ddm_data_freq']}
#   Write variables by module
#       These are automatically enabled when running in batch mode and assigned the same frequency as specified through
#       the Python API, so the user doesn't have to modify any simulation scripts when running a batch
mod_write_vars = {'ViralInfectionVTMModelInputs': ['write_pop_data_freq', 'write_med_diff_data_freq',
                                                   'write_ir_data_freq', 'write_death_data_freq'],
                  'Models.DrugDosingModel.DrugDosingInputs': ['write_ddm_data_freq']}


# todo - test registering unknown input modules with batch workflow
def register_auto_inputs(input_module_name: str, plot_var_names=None, write_var_names=None):
    """
    Registers plot and writing variables with batch run lib for automatic tasks in batch run mode with non-standard
    input modules
    :param input_module_name: string name of input module
    :param plot_var_names: list of string names of plotting variables to automatically disable in batch run mode
    :param write_var_names: list of string names of writing variables to automatically set in batch run mode
    :return: None
    """
    if plot_var_names is not None:
        if input_module_name not in mod_plot_vars.keys():
            mod_plot_vars[input_module_name] = []
        for p in plot_var_names:
            assert isinstance(p, str), f'Specification of variable {p} must be a string'
            if p not in mod_plot_vars[input_module_name]:
                mod_plot_vars[input_module_name].append(p)

    if write_var_names is not None:
        if input_module_name not in mod_write_vars.keys():
            mod_write_vars[input_module_name] = []
        for w in write_var_names:
            assert isinstance(w, str), f'Specification of variable {w} must be a string'
            if w not in mod_write_vars[input_module_name]:
                mod_write_vars[input_module_name].append(w)


def reset_auto_inputs(_input_module_name: str):
    """
    Clears currently registered plot and writing variables
    Can be used to customize pre-defined automatic tasks in distribution
    :param _input_module_name: string name of input module
    :return: None
    """
    if _input_module_name in mod_plot_vars.keys():
        mod_plot_vars[_input_module_name] = []
    if _input_module_name in mod_write_vars.keys():
        mod_write_vars[_input_module_name] = []


def append_auto_inputs(_input_dict: dict) -> None:
    """
    Appends auto inputs to an input dictionary
    :param _input_dict: input dictionary
    :return: None
    """
    assert cc3d_input_key in _input_dict.keys()
    _input_dict[cc3d_input_key][cc3d_auto_key] = [{'input_module_name': k,
                                                   'plot_var_names': mod_plot_vars[k],
                                                   'write_var_names': mod_write_vars[k]
                                                   } for k in mod_plot_vars.keys()]


def apply_external_multipliers(calling_module_str, input_module):
    """
    Applies external input multipliers to model inputs
    Multiplier inputs are stored as a dictionary with key *cc3d_input_key*
    :param calling_module_str: string name of calling module (e.g., 'ViralInfectionVTMSteppables')
    :param input_module: module containing input parameters (e.g., ViralInfectionVTMModelInputs)
    :return:
    """
    from cc3d.CompuCellSetup import persistent_globals
    input_object = persistent_globals.input_object
    if not isinstance(input_object, dict) or cc3d_input_key not in input_object.keys():
        print('No valid batch input found.')
        return

    if cc3d_auto_key in input_object[cc3d_input_key].keys():
        for el in input_object[cc3d_input_key][cc3d_auto_key]:
            reset_auto_inputs(el['input_module_name'])
            register_auto_inputs(**el)
        input_object[cc3d_input_key].pop(cc3d_auto_key)

    for k, v in input_object[cc3d_input_key].items():
        if k != '__param_desc__' and k != cc3d_batch_key:
            try:
                _set_imported_var(calling_module_str, k, v * getattr(input_module, k))
            except AttributeError:
                pass
        elif k == cc3d_batch_key:
            out_freq = v['out_freq']
            if input_module.__name__ in mod_plot_vars.keys():
                for p in mod_plot_vars[input_module.__name__]:
                    _set_imported_var(calling_module_str, p, 0)
            if input_module.__name__ in mod_write_vars.keys():
                for w in mod_write_vars[input_module.__name__]:
                    _set_imported_var(calling_module_str, w, out_freq)


def _set_imported_var(calling_module, var_name, var_val):
    if isinstance(var_name, str):
        return _set_imported_var(sys.modules[calling_module], var_name.split('.'), var_val)
    if len(var_name) == 1:
        setattr(calling_module, var_name[0], var_val)
        return
    return _set_imported_var(getattr(calling_module, var_name[0]), var_name[1:], var_val)


class _MoveDirProcess(multiprocessing.Process):
    def __init__(self, _src_dir, _tgt_dir):
        """
        Process to asynchronously move one directory to another
        :param _src_dir: directory to move
        :param _tgt_dir: location of resulting move
        """
        super().__init__()
        self._src_dir = _src_dir
        self._tgt_dir = _tgt_dir

    def run(self):
        shutil.move(self._src_dir, self._tgt_dir)


def move_dir_async(_src_dir, _tgt_dir) -> None:
    p = _MoveDirProcess(_src_dir, _tgt_dir)
    p.start()
