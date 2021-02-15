"""
Defines module steppables

Steppables
==========

HCVIntegrator
-------------
Description: implements integrated model

Usage: In ViralInfectionVTM.py, add the following

from Models.HCVIntegrated.HCVSteppables import HCVIntegrator

CompuCellSetup.register_steppable(steppable=HCVIntegrator(frequency=1))

HCVCellsInitializerSteppable
----------------------------
Description: changes all initially infected cells' intracellular state to conditions described in Dahari, Harel,
et. al., where infection begins with cytoplasmic viral RNA molecules

Usage: In ViralInfectionVTM.py, add the following

from Models.HCVIntegrated.HCVSteppables import HCVCellsInitializerSteppable

CompuCellSetup.register_steppable(steppable=HCVCellsInitializerSteppable(frequency=1))
"""

from Models.SegoAponte2020 import ViralInfectionVTMLib
from Models.SegoAponte2020.ViralInfectionVTMModelInputs import s_to_mcs
from Models.SegoAponte2020.ViralInfectionVTMSteppables import CellsInitializerSteppable
from Models.SegoAponte2020 import ViralInfectionVTMBasePy

from . import HCVInputs

ViralInfectionVTMSteppableBasePy = ViralInfectionVTMBasePy.ViralInfectionVTMSteppableBasePy

module_prefix = 'ihcv_'

integrator_steppable_key = module_prefix + 'hcv_integrator'  # HCVIntegrator
cell_initializer_key = module_prefix + 'cell_initializer'  # HCVCellsInitializerSteppable


def hcv_viral_replication_model_string(_unpacking_rate, _replicating_rate, _r_half, _translating_rate, _packing_rate,
                                       _secretion_rate, _u_ini, _r_ini, _p_ini, _a_ini, _uptake=0):
    """
    Antimony model string generator for main framework coupled with genomic replication model of hepatitis C virus from
    Dahari, Harel, et al. "Mathematical modeling of subgenomic hepatitis C virus replication in Huh-7 cells."
    Journal of virology 81.2 (2007): 750-760.
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
    model_string = f"""model {ViralInfectionVTMLib.vr_model_name}()
      -> U ; Uptake
    U -> R ; unpacking_rate * U;
    #  -> R ; replicating_rate * r_half * R / (r_half + R);  // Genomic replication from base model; replaced by E1
     -> P ; {HCVInputs.translating_rate} * R;
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

    // Integrated HCV

    // Conversion factors
    s_n = {HCVInputs.virus_from_ul} ;
    s_t = {s_to_mcs / 60 / 60} ;

    // Equations
    E1: -> R ; (k2*TC + kpout*RP)/s_n - k1*RIBO*R - kpin*R - upcyt*R ;
    E2: -> TC    ; k1*RIBO*(R*s_n) - k2*TC - utc*TC ;
    E3: -> PCYT     ; k2*TC - kc*PCYT ;
    E4: -> ECYT  ; kc*PCYT - kein*ECYT - uecyt*ECYT  ;
    E5: -> RP    ; - k3*RP*E + k4p*RIDS + kpin*(R*s_n) - (kpout+ up)*RP;
    E6: -> RDS   ; k4m*RIP + k4p*RIDS - k5*RDS*E - uds*RDS; 
    E7: -> E     ; kein*ECYT + k4m*RIP + k4p*RIDS - k3*RP*E - k5*RDS*E - ue*E;
    E8: -> RIP   ; k3*RP*E - k4m*RIP - uip*RIP;
    E9: -> RIDS  ; k5*RDS*E - k4p*RIDS - uids*RIDS  ;

    RIBO := RIBOTOT - TC

    // Parameters
    k1 = {HCVInputs.k1} * s_t ;
    k2 = {HCVInputs.k2} * s_t ;
    kc = {HCVInputs.kc} * s_t ;
    kpin = {HCVInputs.kpin} * s_t ;
    kpout = {HCVInputs.kpout} * s_t ;
    kein = {HCVInputs.kein} * s_t ;
    k3 = {HCVInputs.k3} * s_t ;
    k4p = {HCVInputs.k4p} * s_t ;
    k4m = {HCVInputs.k4m} * s_t ;
    k5 = {HCVInputs.k5} * s_t ;
    upcyt = {HCVInputs.upcyt} * s_t ;
    up = {HCVInputs.up} * s_t ;
    uds = {HCVInputs.uds} * s_t ;
    uip = {HCVInputs.uip} * s_t ;
    uids = {HCVInputs.uids} * s_t ;
    utc = {HCVInputs.utc} * s_t ;
    ue = {HCVInputs.ue} * s_t ;
    uecyt = {HCVInputs.uecyt} * s_t ;

    // Initial Conditions
    E = 0.0 ;
    ECYT = 0.0 ;
    PCYT = 0.0 ;
    RDS = 0.0 ;
    RIDS = 0.0 ;
    RIP = 0.0 ;
    RP = 0.0 ;
    TC = 0.0 ;
    RIBOTOT  = {HCVInputs.RIBOTOT} ;
    end"""
    return model_string


ihcv_model_vars = ["E", "ECYT", "PCYT", "RDS", "RIBO", "RIDS", "RIP", "RP", "TC"]


class HCVIntegrator(ViralInfectionVTMSteppableBasePy):
    """Integrates HCV model with module in SegoAponte2020"""

    unique_key = integrator_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.set_viral_replication_model(hcv_viral_replication_model_string)


class HCVCellsInitializerSteppable(CellsInitializerSteppable):
    """Initializes cells with initial HCV model data"""

    unique_key = cell_initializer_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self._init_rna = HCVInputs.init_rna
        self._virus_from_ul = HCVInputs.virus_from_ul

    def start(self):
        """
        Called once to initialize simulation
        """
        super().start()
        for cell in self.cell_list_by_type(self.infected_type_id):
            var_unpacking = ViralInfectionVTMLib.vr_cell_dict_to_sym[ViralInfectionVTMLib.vrm_unpacking]
            getattr(cell.sbml, self.vr_model_name)[var_unpacking] = 0

            var_replicating = ViralInfectionVTMLib.vr_cell_dict_to_sym[ViralInfectionVTMLib.vrm_replicating]
            init_replicating = self.init_rna / self.virus_from_ul
            getattr(cell.sbml, self.vr_model_name)[var_replicating] = init_replicating

    @property
    def init_rna(self):
        """
        Initial number of RNA molecules in initially infected cells
        """
        return self._init_rna

    @init_rna.setter
    def init_rna(self, _val: float):
        if _val <= 0:
            raise ValueError("Value must be non-negative")
        self._init_rna = _val

    @property
    def virus_from_ul(self):
        """
        Conversion from unitless original model quantities to HCV model units
        """
        return self._virus_from_ul

    @virus_from_ul.setter
    def virus_from_ul(self, _val: float):
        if _val <= 0:
            raise ValueError("Value must be positive")
        self._virus_from_ul = _val
