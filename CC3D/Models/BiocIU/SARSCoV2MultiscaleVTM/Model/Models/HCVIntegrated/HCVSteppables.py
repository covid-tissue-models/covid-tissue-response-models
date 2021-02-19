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
from Models.SegoAponte2020.ViralInfectionVTMSteppables import CellsInitializerSteppable
from Models.SegoAponte2020 import ViralInfectionVTMBasePy

from . import HCVInputs

ViralInfectionVTMSteppableBasePy = ViralInfectionVTMBasePy.ViralInfectionVTMSteppableBasePy

module_prefix = 'ihcv_'

integrator_steppable_key = module_prefix + 'hcv_integrator'  # HCVIntegrator
cell_initializer_key = module_prefix + 'cell_initializer'  # HCVCellsInitializerSteppable


ihcv_model_vars = ["E", "ECYT", "PCYT", "RDS", "RIBO", "RIDS", "RIP", "RP", "TC"]


class HCVIntegrator(ViralInfectionVTMSteppableBasePy):
    """Integrates HCV model with module in SegoAponte2020"""

    unique_key = integrator_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.set_viral_replication_model(self.hcv_viral_replication_model_string_gen())

        self._virus_from_ul = HCVInputs.virus_from_ul
        self._hcv_params = {'kc': HCVInputs.kc,
                            'kein': HCVInputs.kein,
                            'kpin': HCVInputs.kpin,
                            'kpout': HCVInputs.kpout,
                            'k1': HCVInputs.k1,
                            'k2': HCVInputs.k2,
                            'k3': HCVInputs.k3,
                            'k4p': HCVInputs.k4p,
                            'k4m': HCVInputs.k4m,
                            'k5': HCVInputs.k5,
                            'ribotot': HCVInputs.RIBOTOT,
                            'up': HCVInputs.up,
                            'upcyt': HCVInputs.upcyt,
                            'uds': HCVInputs.uds,
                            'ue': HCVInputs.ue,
                            'uecyt': HCVInputs.uecyt,
                            'uids': HCVInputs.uids,
                            'uip': HCVInputs.uip,
                            'utc': HCVInputs.utc
                            }

    def hcv_viral_replication_model_string_gen(self):
        """
        Returns antimony model string generator function for hooking into SegoAponte2020 module

        Generator function returns model with parameters according to internal data

        :return: model string generator function
        """
        def hcv_viral_replication_model_string(_unpacking_rate, _replicating_rate, _r_half, _translating_rate,
                                               _packing_rate, _secretion_rate, _u_ini, _r_ini, _p_ini, _a_ini,
                                               _uptake=0):
            """
            Antimony model string generator for main framework coupled with genomic
            replication model of hepatitis C virus from

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
              -> P ; translating_rate * R;
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
            s_n = {self.virus_from_ul} ;
            s_t = {self.step_period / 60 / 60} ;
        
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
            kc = {self._hcv_params['kc']} * s_t ;
            kein = {self._hcv_params['kein']} * s_t ;
            kpin = {self._hcv_params['kpin']} * s_t ;
            kpout = {self._hcv_params['kpout']} * s_t ;
            k1 = {self._hcv_params['k1']} * s_t ;
            k2 = {self._hcv_params['k2']} * s_t ;
            k3 = {self._hcv_params['k3']} * s_t ;
            k4p = {self._hcv_params['k4p']} * s_t ;
            k4m = {self._hcv_params['k4m']} * s_t ;
            k5 = {self._hcv_params['k5']} * s_t ;
            uds = {self._hcv_params['uds']} * s_t ;
            ue = {self._hcv_params['ue']} * s_t ;
            uecyt = {self._hcv_params['uecyt']} * s_t ;
            uids = {self._hcv_params['uids']} * s_t ;
            uip = {self._hcv_params['uip']} * s_t ;
            up = {self._hcv_params['up']} * s_t ;
            upcyt = {self._hcv_params['upcyt']} * s_t ;
            utc = {self._hcv_params['utc']} * s_t ;
        
            // Initial Conditions
            E = 0.0 ;
            ECYT = 0.0 ;
            PCYT = 0.0 ;
            RDS = 0.0 ;
            RIDS = 0.0 ;
            RIP = 0.0 ;
            RP = 0.0 ;
            TC = 0.0 ;
            RIBOTOT  = {self._hcv_params['ribotot']} ;
            end"""
            return model_string
        return hcv_viral_replication_model_string

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

    def set_hcv_param(self, **kwargs):
        """
        Set parameters in hcv model for newly created cells

        Keywords are

        - kc: Viral polyprotein cleavage (1/h)
        - kein: ECYT transport into VMS (1/h)
        - kpout: RP transport into cytoplasm (1/h)
        - kpin: R transport into VMS (1/h)
        - k1: Tc formation (1/h)
        - k2: Nascent NS polyprotein translation (1/h)
        - k3: RIP formation (1/h)
        - k4p: RP synthesis (1/h)
        - k4m: RDS synthesis (1/h)
        - k5: RIDS formation (1/h)
        - ribotot: Total ribosomes
        - upcyt: R degradation (1/h)
        - uds: RDS degradation (1/h)
        - ue: E degradation (1/h)
        - uecyt: ECYT degradation (1/h)
        - uids: RIDS degradation (1/h)
        - uip: RIP degradation (1/h)
        - up: RP degradation (1/h)
        - utc: TC degradation (1/h)

        Any model parameter can be set, so long as the value is non-negative

        :param kwargs: keyword argument values
        :return: None
        """
        for k, v in kwargs.items():
            if k not in self._hcv_params.keys():
                raise AttributeError(f'Unrecognized parameter ({k} = {v})')
            elif v < 0.0:
                raise ValueError(f'Value must be non-negative ({k} = {v})')
            self._hcv_params[k] = v


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
