# This is a library for the coronavirus viral infection modeling project using CompuCell3D
# by the Biocomplexity Institute at Indiana University

import os
import sys

from cc3d.cpp import CompuCell
import numpy as np

# Key to mcs value when a cell was created
new_cell_mcs_key = 'new_cell_mcs'

# Key to cell dictionary boolean for whether an instance of the viral replication model has been loaded
vrl_key = 'viral_replication_loaded'

# Name of Antimony/SBML model of viral replication
vr_model_name = 'viralReplication'

# Mapping from CellG instance dictionary keys to Antimony/SBML viral replication model symbols
vr_cell_dict_to_sym = {'Unpacking': 'U',
                       'Replicating': 'R',
                       'Packing': 'P',
                       'Assembled': 'A'}

# Key to cell dictionary boolean for whether an instance of the viral internalization model has been loaded
vil_key = 'viral_internalization_loaded'

# Name of Antimony/SBML viral internalization model
vi_model_name = 'viralInternalization'

# Mapping from CellG instance dictionary keys to Antimony/SBML viral internalization model symbols
vi_cell_dict_to_sym = {'Unbound_Receptors': 'R',
                       'Surface_Complexes': 'VR',
                       'Internalized_Complexes': 'Vi'}


# Name of Antimony/SBML model of immune cell recruitment
ir_model_name = 'immuneRecruitment'

# Key to reference of ImmuneRecruitmentSteppable instance in shared global dictionary
ir_steppable_key = 'ir_steppable'

# Key to reference of SimDataSteppable instance in shared global dictionary
simdata_steppable_key = 'simdata_steppable'


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


def viral_internalization_model_string(_kon, _koff, _intern_rate, _ve_ini=0, _r_ini=0, _vr_ini=0, _vi_ini=0, _ve_src=0):
    """
    dVe/dt = -kon*Ve*R + koff*VR + VeSrc
    dR/dt = -kon*Ve*R + koff*VR
    dVR/dt = kon*Ve*R - koff*VR - intern_rate*VR
    dVi/dt = intern_rate*VR
    Derived by J. Aponte-Serrano and J. Mathur
    :param _kon: association rate constant of extracellular virus particles and unbound cell receptors
    :param _koff: dissasociation rate constant of virus-receptor surface complex
    :param _intern_rate: internalization rate of virus-receptor surface complex
    :param _ve_ini: initial number of extracellular virus particles
    :param _r_ini: initial number of unbound cell receptors
    :param _vr_ini: initial number of virus-receptor surface complexes
    :param _vi_ini: initial number of internalized virus particles
    :param _ve_src: incoming extracellular virus particles from viral field
    :return: None
    """
    model_string = """model {}()
          -> Ve; VeSrc;
        Ve + R  -> VR ; kon * Ve * R ;
        VR -> Ve + R  ; koff * VR ;
        VR -> Vi ; intern_rate * VR ;
        kon = {};
        koff = {};
        intern_rate = {};
        VeSrc = {};
        Ve = {};
        R = {};
        VR = {};
        Vi = {};
        end""".format(vi_model_name, _kon, _koff, _intern_rate, _ve_src, _ve_ini, _r_ini, _vr_ini, _vi_ini)
    return model_string


def immune_recruitment_model_string(_add_rate, _sub_rate, _delay_rate, _decay_rate, _total_ck=0, _num_imm=0, _s_ini=0):
    """
    dS/dt = addRate - subRate * numImmuneCells + delayRate * totalCytokine - decayRate * S
    The probability of adding an immume cell is non-zero for S > 0
    The probabiilty of removing an immune cell is non-zero for S < 0
    Derived in part thanks to J. Toledo
    :param _add_rate: addition rate
    :param _sub_rate: substraction rate
    :param _delay_rate: delay rate
    :param _decay_rate: decay rate
    :param _total_ck: total cytokine signal
    :param _num_imm: total number of immune cells
    :param _s_ini: initial value of state variable *S*
    :return: None
    """
    model_string = """model {}()
          -> S ; addRate + delayRate * totalCytokine;
        S ->   ; subRate * numImmuneCells + decayRate * S;
        addRate = {};
        subRate = {};
        delayRate = {};
        decayRate = {};
        numImmuneCells = {};
        totalCytokine = {};
        S = {};
        end""".format(ir_model_name, _add_rate, _sub_rate, _delay_rate, _decay_rate, _total_ck, _num_imm, _s_ini)
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
    Loads state variables from viral replication model SBML into cell dictionary
    :param cell: cell for which to load state variables from SBML into cell dictionary
    :return: None
    """
    assert cell.dict[vrl_key]
    for k, v in vr_cell_dict_to_sym.items():
        cell.dict[k] = getattr(cell.sbml, vr_model_name)[v]


def pack_viral_internalization_variables(cell):
    """
    Loads state variables from viral internalization model SBML into cell dictionary
    :param cell: cell for which to load state variables from SBML into cell dictionary
    :return: None
    """
    assert cell.dict[vil_key]
    for k, v in vi_cell_dict_to_sym.items():
        cell.dict[k] = getattr(cell.sbml, vi_model_name)[v]


def reset_viral_replication_variables(cell):
    """
    Sets state variables from viral replication model in cell dictionary to zero
    :param cell: cell for which to set state variables in cell dictionary to zero
    :return: None
    """
    cell.dict['Uptake'] = 0
    cell.dict['Secretion'] = 0
    for k in vr_cell_dict_to_sym.keys():
        cell.dict[k] = 0


def reset_viral_internalization_variables(cell):
    """
    Sets state variables from viral internalization model in cell dictionary to zero
    :param cell: cell for which to set state variables in cell dictionary to zero
    :return: None
    """
    for k in vi_cell_dict_to_sym.keys():
        cell.dict[k] = 0


def internalize_viral_particles(cell, vi_step_size):
    """
    Moves internalized viral particles from internalization model to replication model for a cell
    :param cell: cell for which to perform transfer of internalized viral particles
    :param vi_step_size: time step size of viral internalization model
    :return: None
    """
    assert cell.dict[vrl_key] and cell.dict[vil_key]
    vi_sbml = getattr(cell.sbml, vi_model_name)
    intern_vir = vi_sbml['Vi']
    vi_sbml['Vi'] = 0.0
    set_viral_replication_cell_uptake(cell, intern_vir / vi_step_size)


def step_sbml_viral_internalization_cell(cell, vi_step_size, ve_tot=0):
    """
    Steps viral internalization SBML model for a cell
    :param cell: cell with a SBML model to step
    :param vi_step_size: time step size of viral internalization model
    :param ve_tot: total environmental viral amount to pass to SBML model
    :return: extracellular viral amount after integration of SBML model
    """
    assert cell.dict[vil_key]
    vi_sbml = getattr(cell.sbml, vi_model_name)
    vi_sbml['VeSrc'] = ve_tot / vi_step_size
    step_sbml_model_cell(cell, vi_model_name)
    extern_vir = vi_sbml['Ve']
    vi_sbml['Ve'] = 0.0
    return extern_vir


def get_assembled_viral_load_inside_cell(cell, sbml_rate):
    return sbml_rate*cell.dict['Uptake'] + cell.dict['Assembled']
