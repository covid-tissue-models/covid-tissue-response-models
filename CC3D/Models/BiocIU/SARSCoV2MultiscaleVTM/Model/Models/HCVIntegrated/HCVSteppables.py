# Integrated model of hepatitis c virus
# Written by T.J. Sego, Ph.D.
# Genomic replication from the main framework is replaced by a model of hepatitis c virus replication from
#
#   Dahari, Harel, et al. "Mathematical modeling of subgenomic hepatitis C virus replication in Huh-7 cells."
#   Journal of virology 81.2 (2007): 750-760.
#
# Model parameters of the hepatitis c virus model are specified in HCVInputs.py
#
# HCVIntegrator
#   Description: implements integrated model
#   Usage:
#       In ViralInfectionVTM.py, add the following
#           from Models.HCVIntegrated.HCVSteppables import HCVIntegrator
#           CompuCellSetup.register_steppable(steppable=HCVIntegrator(frequency=1))
# HCVDataSteppable
#   Description: performs data tracking
#   Usage:
#       In ViralInfectionVTM.py, add the following
#           from Models.HCVIntegrated.HCVSteppables import HCVDataSteppable
#           CompuCellSetup.register_steppable(steppable=HCVDataSteppable(frequency=1))

import os
import sys

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
import ViralInfectionVTMLib
from ViralInfectionVTMModelInputs import s_to_mcs
from ViralInfectionVTMSteppables import SimDataSteppable
from ViralInfectionVTMSteppableBasePy import ViralInfectionVTMSteppableBasePy

from . import HCVInputs


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
    def __init__(self, frequency=1):
        super().__init__(frequency)
        self.set_viral_replication_model(hcv_viral_replication_model_string)


class HCVDataSteppable(ViralInfectionVTMSteppableBasePy):
    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self.ihcv_data_win = None
        self.ihcv_data_path = None
        self.ihcv_data = dict()

        self.plot_ihcv_data = HCVInputs.plot_ihcv_data_freq > 0
        self.write_ihcv_data = HCVInputs.write_ihcv_data_freq > 0

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

    def start(self):
        if self.plot_ihcv_data:
            self.ihcv_data_win = self.add_new_plot_window(title='Integrated HCV',
                                                          x_axis_title='MonteCarlo Step (MCS)',
                                                          y_axis_title='Variables', x_scale_type='linear',
                                                          y_scale_type='linear',
                                                          grid=False,
                                                          config_options={'legend': True})

            colors = ['blue', 'red', 'green', 'yellow', 'white', 'cyan', 'purple', 'magenta', 'darkgreen']
            for x in range(len(ihcv_model_vars)):
                self.ihcv_data_win.add_plot(ihcv_model_vars[x], style='Dots', color=colors[x], size=5)

        if self.write_ihcv_data:
            from pathlib import Path
            self.ihcv_data_path = Path(self.output_dir).joinpath('ihcv_data.dat')
            with open(self.ihcv_data_path, 'w'):
                pass

    def step(self, mcs):
        if self.simdata_steppable is None:
            self.simdata_steppable = self.shared_steppable_vars[ViralInfectionVTMLib.simdata_steppable_key]

        plot_ihcv_data = self.plot_ihcv_data and mcs % HCVInputs.plot_ihcv_data_freq == 0
        write_ihcv_data = self.write_ihcv_data and mcs % HCVInputs.write_ihcv_data_freq == 0

        vrm_tracked_cell = self.simdata_steppable.vrm_tracked_cell
        if vrm_tracked_cell is not None and (plot_ihcv_data or write_ihcv_data):
            cell_sbml = getattr(vrm_tracked_cell.sbml, self.vr_model_name)
            if plot_ihcv_data:
                [self.ihcv_data_win.add_data_point(x, mcs, cell_sbml[x]) for x in ihcv_model_vars]

            if write_ihcv_data:
                self.ihcv_data[mcs] = [vrm_tracked_cell.id]
                [self.ihcv_data[mcs].append(cell_sbml[x] for x in ihcv_model_vars)]

        # Flush outputs at quarter simulation lengths
        if mcs >= int(self.simulator.getNumSteps() / 4 * self.__flush_counter):
            self.flush_stored_outputs()
            self.__flush_counter += 1

    def on_stop(self):
        self.finish()

    def finish(self):
        self.flush_stored_outputs()

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage
        :return: None
        """
        if self.write_ihcv_data:
            with open(self.ihcv_data_path, 'a') as fout:
                fout.write(SimDataSteppable.data_output_string(self.ihcv_data))
                self.ihcv_data.clear()
