# This is a library for the viral infection modeling project using CompuCell3D
# by the Biocomplexity Institute at Indiana University

# Module prefix
from . import module_prefix

# Key to mcs value when a cell was created
new_cell_mcs_key = module_prefix + 'new_cell_mcs'

# Key to cell dictionary boolean for whether an instance of the viral replication model has been loaded
vrl_key = module_prefix + 'viral_replication_loaded'

# Name of Antimony/SBML model of viral replication
vr_model_name = module_prefix + 'viralReplication'

# Mapping from CellG instance dictionary keys to Antimony/SBML viral replication model symbols
vrm_uptake = module_prefix + 'Uptake'
vrm_unpacking = module_prefix + 'Unpacking'
vrm_replicating = module_prefix + 'Replicating'
vrm_packing = module_prefix + 'Packing'
vrm_assembled = module_prefix + 'Assembled'
vrm_secretion = module_prefix + 'Secretion'
vr_cell_dict_to_sym = {vrm_unpacking: 'U',
                       vrm_replicating: 'R',
                       vrm_packing: 'P',
                       vrm_assembled: 'A'}

# Name of Antimony/SBML model of immune cell recruitment
ir_model_name = module_prefix + 'immuneRecruitment'

# Unique keys for module steppables
cell_initializer_steppable_key = module_prefix + 'cell_initializer_key'  # CellsInitializerSteppable
virus_field_initializer_key = module_prefix + 'virus_field_initializer_steppable'  # VirusFieldInitializerSteppable
vrm_steppable_key = module_prefix + 'vrm_steppable'  # ViralReplicationSteppable
vim_steppable_key = module_prefix + 'vim_steppable'  # ViralInternalizationSteppable
vrs_steppable_key = module_prefix + 'vrs_steppable'  # ViralSecretionSteppable
immune_killing_steppable_key = module_prefix + 'immune_killing_steppable'  # ImmuneCellKillingSteppable
chemotaxis_steppable_key = module_prefix + 'chemotaxis_steppable'  # ChemotaxisSteppable
immune_seeding_steppable_key = module_prefix + 'immune_seeding_steppable'  # ImmuneCellSeedingSteppable
simdata_steppable_key = module_prefix + 'simdata_steppable'  # SimDataSteppable
cytokine_secretion_steppable_key = module_prefix + 'cytokine_secretion_steppable'  # CytokineProductionAbsorptionSteppable
ir_steppable_key = module_prefix + 'ir_steppable'  # ImmuneRecruitmentSteppable
oxidation_steppable_key = module_prefix + 'oxidation_steppable'  # oxidationAgentModelSteppable

# Model parameter keys in cell dictionary
unbound_receptors_cellg_key = module_prefix + 'Receptors'
ck_production_cellg_key = module_prefix + 'ck_production'
ck_consumption_cellg_key = module_prefix + 'ck_consumption'
tot_ck_upt_cellg_key = module_prefix + 'tot_ck_upt'
activated_cellg_key = module_prefix + 'activated'
time_activation_cellg_key = module_prefix + 'time_activation'
oxi_killed_cellg_key = module_prefix + 'oxi_killed'


# todo: Generalize Antimony model string generator for general use
def viral_replication_model_string(_unpacking_rate, _replicating_rate, _r_half, _translating_rate, _packing_rate,
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
    :param _r_half: Value of R at which the replication rate is half max
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
    model_string = f"""model {vr_model_name}()
      -> U ; Uptake
    U -> R ; unpacking_rate * U;
      -> R ; replicating_rate * r_half * R / (r_half + R);
    R -> P ; translating_rate * R;
    P -> A ; packing_rate * P;
    A -> Secretion ; secretion_rate * A;

    unpacking_rate = {_unpacking_rate};
    replicating_rate = {_replicating_rate};
    r_half = {_r_half};
    translating_rate = {_translating_rate};
    packing_rate = {_packing_rate};
    secretion_rate = {_secretion_rate};
    U = {_u_ini};
    R = {_r_ini};
    P = {_p_ini};
    A = {_a_ini};
    Uptake = {_uptake};
    Secretion = 0;
    end"""
    return model_string


def immune_recruitment_model_string(_add_rate, _sub_rate, _delay_rate, _decay_rate, _total_ck=0, _num_imm=0, _s_ini=0):
    """
    dS/dt = addRate - subRate * numImmuneCells + totalCytokine / delayRate - decayRate * S

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
    model_string = f"""model {ir_model_name}()
          -> S ; addRate + totalCytokine / delayRate;
        S ->   ; subRate * numImmuneCells + decayRate * S;
        addRate = {_add_rate};
        subRate = {_sub_rate};
        delayRate = {_delay_rate};
        decayRate = {_decay_rate};
        numImmuneCells = {_num_imm};
        totalCytokine = {_total_ck};
        S = {_s_ini};
        end"""
    return model_string


def step_sbml_model_cell(cell, sbml_model_name=vr_model_name):
    """
    Steps SBML model for a cell

    :param cell: cell with a SBML model to step
    :param sbml_model_name: name of SBML model to step
    :return: None
    """
    from cc3d.cpp import CompuCell
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


def reset_viral_replication_variables(cell):
    """
    Sets state variables from viral replication model in cell dictionary to zero

    :param cell: cell for which to set state variables in cell dictionary to zero
    :return: None
    """
    cell.dict[vrm_uptake] = 0
    cell.dict[vrm_secretion] = 0
    for k in vr_cell_dict_to_sym.keys():
        cell.dict[k] = 0


def get_assembled_viral_load_inside_cell(cell, sbml_rate):
    return sbml_rate*cell.dict[vrm_uptake] + cell.dict[vrm_assembled]
