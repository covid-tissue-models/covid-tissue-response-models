# -*- coding: utf-8 -*-
"""
Created on Wed Sep 22 16:03:29 2021

@author: Juliano
"""
import os
import numpy as np
import pandas as pd
import json
from scipy.signal import argrelextrema, find_peaks

def get_viral_load_column(path, usecols=[1]):
    with open(path) as csvf:
        df = pd.read_csv(path, header=None, usecols=usecols)

    return df.to_numpy()

def get_sets_of_param_value(parameter_name, value, inputs):
    sets = []
    for idx, in_dict in enumerate(inputs):
        this_value = in_dict['__input_dict__'][parameter_name]
        if this_value == value:
            sets.append(f'set_{idx}')
    return sets

def get_parameters_by_name(names, indict):
    parameters = {}
    for key, value in indict.items():
        if key in names:
            parameters[key] = value

    return parameters

def get_sim_inputs(path):
    with open(os.path.join(path, 'batch_status.json'), 'r') as batch_json:
        sim_inputs = json.load(batch_json)['sim_input']
    return sim_inputs

def get_all_runs_data(file_name, path, runs, usecols=[1]):
    
    tcs = np.zeros((8001, len(runs)))
    
    for i, run in enumerate(runs):
        if not os.path.isfile(os.path.join(path,run,file_name)):
            file_name = file_name.split('.')[0]+'.dat'
        tc =  get_viral_load_column(os.path.join(path,run,file_name), usecols=usecols).flatten()
        if len(tc) < 8001:
            print('SHORT FILE', run)
            tcs[:len(tc),i] = tc[:]
        else:
            tcs[:,i] = get_viral_load_column(os.path.join(path,run,file_name), usecols=usecols).flatten()
    return tcs

def get_slope(vir_median, fast_containment_thrs):
    log_vir = np.log10(vir_median[fast_containment_thrs:])
    minima = argrelextrema(log_vir, np.less)
    maxima, _ = find_peaks(log_vir)
    
    if len(minima[0]) and len(maxima):
        slope = log_vir[maxima[-1]] - log_vir[minima[0][0]] 
    elif len(maxima):
        slope = log_vir[maxima[-1]] - log_vir[fast_containment_thrs] 
    elif len(minima):
        slope = log_vir[-1] - log_vir[minima[0][0]] 
    else:
        slope = log_vir[-1] - log_vir[fast_containment_thrs]
    return slope

def peaks_slope(vir_median, median_time_start_idx, fast_containment_thrs):
    log_vir = np.log10(vir_median[median_time_start_idx:])
    maxima_idx, _ = find_peaks(log_vir)
    if len(maxima_idx):
        slope, _ =  np.polyfit(maxima_idx, log_vir[maxima_idx], 1)
        return slope
    else:
        
        slope = log_vir[-1] - np.log10(np.min(vir_median[:fast_containment_thrs//2 - median_time_start_idx]))
        return slope


def determine_tendency(slope, chronic_thr, auc_from_14, auc_thr):
    if slope < -chronic_thr:
        return 'containment'
    elif abs(slope) < chronic_thr or auc_from_14 < auc_thr:
        return 'chronic'
    elif slope < 0:
        return 'containment'
    else:
        return 'runaway'

def get_median_treatment_start(file, path, runs, first_dose, init_infected=5, time_conv= 5 * 60 / 60  / 60  / 24):
    inf_tcs = get_all_runs_data(file, path, runs, usecols=[2])
    infections_tcs = np.diff(inf_tcs, axis=0) 
    new_inf =  np.where(infections_tcs>=0, infections_tcs, 0)
    ever_inf = np.zeros(new_inf.shape)
    ever_inf[0,:] = init_infected + new_inf[0,:]
    for i, row in enumerate(new_inf[1:]):
        ever_inf[i+1,:] = ever_inf[i,:] + row[:]
    treat_timer_start = np.where(ever_inf>=10, True, False)
    time_start_idx = []
    for col in treat_timer_start.T:
        time_start_idx.append(np.argmax(col) + first_dose/time_conv)
    median_time_start_idx = np.median(time_start_idx)
    median_time_start = median_time_start_idx*time_conv
    median_time_start_idx = round(median_time_start_idx)
    return median_time_start, median_time_start_idx

def determine_outcome(s, runs, path_of_batches, file_vir, file_pop, file_auc, time_conv, slope_thr, first_dose,
                      thrs_vir=1.3, chronic_thr = 1e-4):
    
    healthy_tcs = get_all_runs_data(file_pop, os.path.join(path_of_batches,s),runs, usecols=[1])
    healthy_median = np.median(healthy_tcs, axis=1)
    median_time_start, median_time_start_idx = get_median_treatment_start(
        file_pop, os.path.join(path_of_batches,s),runs, first_dose)
    if healthy_median[-1] <= 10:
        cleared_without_reapearance = cleared_fast = cleared_vir = False
        tendency = 'runaway'
        return cleared_vir, cleared_without_reapearance, cleared_fast, tendency, median_time_start
    
    healthy_at_treat_start = healthy_median[median_time_start_idx]
    
    if healthy_median[-1] < 0.5*healthy_at_treat_start:
        cleared_without_reapearance = cleared_fast = cleared_vir = False
        tendency = 'runaway'
        return cleared_vir, cleared_without_reapearance, cleared_fast, tendency, median_time_start
    fast_containment_thrs = median_time_start_idx + int(14/time_conv)
    
    vir_tcs = get_all_runs_data(file_vir,os.path.join(path_of_batches,s),runs)
        
    vir_median = np.median(vir_tcs, axis=1)
    low_vir_idx = np.where(vir_median[median_time_start_idx:]<thrs_vir)[0]
    cleared_vir = bool(len(low_vir_idx))
    
    if not cleared_vir:
        
        auc_tcs = get_all_runs_data(file_auc,os.path.join(path_of_batches,s),runs)
        
        auc_median = np.median(auc_tcs, axis=1)
        
        auc_from_14 = auc_median[-1] - auc_median[fast_containment_thrs]
        
        cleared_without_reapearance = cleared_fast = cleared_vir
        thrs_vir_secodnary = 1.1*thrs_vir
        not_so_low_vir = np.where(vir_median[median_time_start_idx:]<thrs_vir_secodnary)[0]
        if len(not_so_low_vir):
            return cleared_vir, cleared_without_reapearance, cleared_fast, 'containment', median_time_start
        slope = peaks_slope(vir_median, median_time_start_idx, fast_containment_thrs)
        auc_thr = 300
        tendency = determine_tendency(slope, chronic_thr, auc_from_14, auc_thr)
        
        
        return cleared_vir, cleared_without_reapearance, cleared_fast, tendency, median_time_start
    else:
        high_vir_idx = np.where(vir_median[median_time_start_idx:]>thrs_vir)[0]
        low_high_vir_idx = np.where(vir_median[median_time_start_idx:]<1.1*thrs_vir)[0]
        consecutive_list = np.array(list(range(np.min(low_high_vir_idx), np.max(low_high_vir_idx)+1)))
        cleared_without_reapearance = ((low_high_vir_idx[-1] == (len(vir_median[median_time_start_idx:])-1)) and 
                                           (len(consecutive_list) == len(low_high_vir_idx)))
        cleared_fast = low_vir_idx[0] < fast_containment_thrs and high_vir_idx[-1] < fast_containment_thrs
        # if not cleared_without_reapearance or not cleared_fast:
        if not cleared_fast:
            slope = get_slope(vir_median, fast_containment_thrs)
            tendency = determine_tendency(slope, chronic_thr, 15, 10)
            return cleared_vir, cleared_without_reapearance, cleared_fast, tendency, median_time_start
        else:
            tendency = 'containment'
            return cleared_vir, cleared_without_reapearance, cleared_fast, tendency, median_time_start
        
    
    return None, None, None, None, None

def pick_color(check_list, colors):
    if check_list[0] and check_list[2]:
        return colors[0]
    elif check_list[0] and not check_list[2]:
        return colors[1]
    elif check_list[-1] == 'containment' or (check_list[1] and not check_list[2]):
        return colors[1]
    elif check_list[-1] == 'runaway':
        return colors[-1]
    elif check_list[-1] == 'chronic':
        return colors[-2]
    return None