# This is a library for the coronavirus viral infection modeling project using CompuCell3D
# by the Biocomplexity Institute at Indiana University

import os
import sys

from cc3d.cpp import CompuCell
import numpy as np

# Key to cell dictionary boolean for whether an instance of the viral replication model has been loaded
vrl_key = 'viral_replication_loaded'

# Name of Antimony/SBML model
vr_model_name = 'viralReplication'

# Mapping from CellG instance dictionary keys to Antimony/SBML symbols
vr_cell_dict_to_sym = {'Unpacking': 'U',
                       'Replicating': 'R',
                       'Packing': 'P',
                       'Assembled': 'A'}


# todo: Generalize Antimony model string generator for general use
def viral_replication_model_string(_unpacking_rate, _replicating_rate, _translating_rate, _packing_rate,
                                   _secretion_rate, _u_ini, _r_ini, _p_ini, _a_ini, _uptake=0):
    """
    Antimony model string generator for this project
    To change models, modify according to this structure
    Modifications should be reflected in
      1. the items of the dictionary "vr_cell_dict_to_sym"
      2. the signature of the function "load_viral_replication_model"
    Variable "Uptake" is the uptake variable of the model, and should not be modified
    Variable "Secretion" is the secretion variable of the model, and should not be modified
    :param _unpacking_rate: model unpacking rate
    :param _replicating_rate: model replicating rate
    :param _translating_rate: model translating rate
    :param _packing_rate: model packing rate
    :param _secretion_rate: model secretion rate
    :param _u_ini: initial model U
    :param _r_ini: initial model R
    :param _p_ini: initial model P
    :param _a_ini: initial model A
    :param _uptake: model Uptake
    :return: Antimony model string
    """
    model_string = """model {}()
      -> U ; Uptake
    U -> R ; unpacking_rate * U;
      -> R ; replicating_rate * R / (0.1 + R);
    R -> P ; translating_rate * R;
    P -> A ; packing_rate * P;
    A -> Secretion ; secretion_rate * A;

    unpacking_rate = {};
    replicating_rate = {};
    translating_rate = {};
    packing_rate = {};
    secretion_rate = {};
    U = {};
    R = {};
    P = {};
    A = {};
    Uptake = {};
    Secretion = 0;
    end""".format(vr_model_name,
                  _unpacking_rate, _replicating_rate, _translating_rate, _packing_rate, _secretion_rate,
                  _u_ini, _r_ini, _p_ini, _a_ini, _uptake)
    return model_string


def step_sbml_model_cell(cell, sbml_model_name=vr_model_name):
    """
    Steps SBML model for a cell
    :param cell: cell with a SBML model to step
    :param sbml_model_name: name of SBML model to step
    :return: None
    """
    dict_attrib = CompuCell.getPyAttrib(cell)
    assert 'SBMLSolver' in dict_attrib
    dict_attrib['SBMLSolver'][sbml_model_name].timestep()


def enable_viral_secretion(cell, secretion_rate, _enable: bool = True):
    """
    Enable/disable secretion in state model for a cell
    :param cell: cell for which to enable/disable secretion in state model
    :param secretion_rate: value for which to get the value of state variable "Secretion"
    :param _enable: enables secretion with value *secretion_rate* when True; sets value to zero when False
    :return:
    """
    if _enable:
        getattr(cell.sbml, vr_model_name)['secretion_rate'] = secretion_rate
    else:
        getattr(cell.sbml, vr_model_name)['secretion_rate'] = 0.0


def set_viral_replication_cell_uptake(cell, uptake):
    """
    Sets the current state variable "Uptake" for a cell
    :param cell: cell for which to set the value of state variable "Uptake"
    :param uptake: value to set for state variable "Uptake"
    :return: None
    """
    assert cell.dict[vrl_key]
    getattr(cell.sbml, vr_model_name)['Uptake'] = uptake


def get_viral_replication_cell_secretion(cell):
    """
    Gets the current state variable "Secretion" for a cell
    :param cell: cell for which to get the value of state variable "Secretion"
    :return: value of state variable "Secretion"
    """
    assert cell.dict[vrl_key]
    secr = getattr(cell.sbml, vr_model_name)['Secretion']
    getattr(cell.sbml, vr_model_name)['Secretion'] = 0.0
    return secr


def pack_viral_replication_variables(cell):
    """
    Loads state variables from SBML into cell dictionary
    :param cell: cell for which to load state variables from SBML into cell dictionary
    :return: None
    """
    assert cell.dict[vrl_key]
    for k, v in vr_cell_dict_to_sym.items():
        cell.dict[k] = getattr(cell.sbml, vr_model_name)[v]


def reset_viral_replication_variables(cell):
    """
    Sets state variables in cell dictionary to zero
    :param cell: cell for which to set state variables in cell dictionary to zero
    :return: None
    """
    cell.dict['Uptake'] = 0
    cell.dict['Secretion'] = 0
    for k in vr_cell_dict_to_sym.keys():
        cell.dict[k] = 0
