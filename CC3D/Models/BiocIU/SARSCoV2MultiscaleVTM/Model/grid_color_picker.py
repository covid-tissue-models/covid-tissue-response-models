# -*- coding: utf-8 -*-
"""
Created on Thu Aug  5 14:57:57 2021

@author: Juliano
"""

import os
import numpy as np

from grid_color_picker_functions import get_sim_inputs, get_parameters_by_name, get_sets_of_param_value, \
                                        determine_outcome, get_all_runs_data, pick_color


import json


if __name__ == '__main__':
    
    path_of_sets = r'C:\my_results' # This should point to the folder that contains the set_* folders

    inputs = get_sim_inputs(path_of_sets)
    
    
    first_dose = get_parameters_by_name('first_dose', inputs[0]['__input_dict__'])['first_dose']
    
    sets_to_do = get_sets_of_param_value('first_dose', first_dose, inputs)
    
    runs = ['run_0', 'run_1', 'run_2', 'run_3', 'run_4', 'run_5', 'run_6', 'run_7']
    
    file_vir = 'med_diff_data.csv'
    file_pop = 'pop_data.csv'
    file_auc = 'ddm_total_viral_production_data.csv'
    thrs_vir = 1.3
    
    chronic_thr = 1e-4
    init_infected = 5
    time_conv = 5 * 60 / 60  / 60  / 24
    
    log_slopes = []
    
    
    analysis_dict = {}
    
    
    colors =['#008000', '#0000FF', '#000000', '#FF0000']
    
    names = ['Rapid clearance', 'Slow clearance',
             'Persistent infection', 'Runaway virus']
    
    analysis_dict['comment'] = [['cleared virus', 'cleared without reapearance (unused)', 'fast clearance', 
                                 'tendency (chronic, runaway, containment)'],
                                'color',
                                'vir thrs'
                                'median time start',
                                'median time start + 14']
    
    for s in sets_to_do:
        print(s)
        
        cleared_vir, cleared_without_reapearance, \
            cleared_fast, tendency, median_start_time = determine_outcome(s, runs, path_of_sets, file_vir,
                                                       file_pop, file_auc, time_conv, chronic_thr, first_dose)
        cleared_fast = bool(cleared_fast)
        vir_tcs = get_all_runs_data(file_vir,os.path.join(path_of_sets,s),runs)
        
        vir_median = np.median(vir_tcs, axis=1)
        
        analysis_dict[s] = [[bool(cleared_vir), bool(cleared_without_reapearance), 
                             bool(cleared_fast), tendency],
                            None,
                            thrs_vir,
                            median_start_time,
                            median_start_time+14]
        analysis_dict[s][1] = pick_color(analysis_dict[s][0], colors)
        
    with open(os.path.join(path_of_sets, 'set_colors.json'), 'w+') as f:
        json.dump(analysis_dict, f, indent=4)