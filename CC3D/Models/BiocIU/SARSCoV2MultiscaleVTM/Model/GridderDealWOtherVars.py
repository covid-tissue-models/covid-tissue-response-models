# Written by Juliano F Gianlupi, MS

import os
import json


def get_parameter_names(sim_inputs):
    return list(sim_inputs[0]['__input_dict__'].keys())


def get_folder_parameters(parameter_names, parameters_to_grid):
    return [p for p in parameter_names if p not in parameters_to_grid and '__' not in p]

def get_non_grid_folder_for_this_param_value(parameters_to_grid, value, path):
    return
    folder = get_parameter_names(get_sim_inputs(path))

    # folder = [p for p in folder if value in ]


def get_sim_inputs(path):
    with open(os.path.join(path, 'batch_status.json'), 'r') as batch_json:
        sim_inputs = json.load(batch_json)['sim_input']
    return sim_inputs


def get_sets_of_param_value(parameter_name, value, inputs):
    sets = []
    for idx, in_dict in enumerate(inputs):
        this_value = in_dict['__input_dict__'][parameter_name]
        if this_value == value:
            sets.append(f'set_{idx}')
    return sets


if __name__ == '__main__':
    dir_spatial_results_root_abs = r'D:\Google Drive IU\phdStuff\covid 19 project\ddm results\new ' \
                                   r'PK\correct_ic50\ddm_batch_7'
    inputs = get_sim_inputs(dir_spatial_results_root_abs)

    sets = get_sets_of_param_value('first_dose', 0, inputs)

    a = 10
