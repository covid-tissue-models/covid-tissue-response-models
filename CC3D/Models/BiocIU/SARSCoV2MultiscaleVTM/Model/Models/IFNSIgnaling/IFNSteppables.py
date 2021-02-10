# IFN intracellular signaling and coupled viral replication models
# Written by Josua Aponte-Serrano
# Adds IFN signaling model to cells in the main framework. Replaces viral replication model in the main framework
# Adopted from:
#
# Multiscale Model of RNA Virus Replication and Interferon Responses Reveals Factors Controlling Plaque
# Growth Dynamics" Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego, James A. Glazier and
# Jason E. Shoemaker
#
# Model parameters are specified in IFNInputs.py

import os
import sys

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
import ViralInfectionVTMLib
from ViralInfectionVTMModelInputs import s_to_mcs
from ViralInfectionVTMSteppableBasePy import ViralInfectionVTMSteppableBasePy
from ViralInfectionVTMSteppables import SimDataSteppable
from . import IFNInputs

def IFN_model_string(_H=0.0,_IFNe=0.0,_V=0.0):
    """
    Antimony model string generator for IFN intracellular signaling model adopted from "Multiscale Model
    of RNA Virus Replication and Interferon Responses Reveals Factors Controlling Plaque Growth Dynamics"
    Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego, James A. Glazier and Jason E. Shoemaker
    :return: Antimony model string
    """
    model_string = f'''IFN_model()
        //Equations
        E2a: -> IFN         ; H*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)    ;
        E2b: IFN ->         ; k21*IFN                                           ;
        E4a: -> STATP       ; H*k31*IFNe/(k32+k33*IFNe)                         ;
        E4b: STATP ->       ; t3*STATP                                          ;
        E5a: -> IRF7        ; H*(k41*STATP+k42*IRF7P)                           ;
        E5b: IRF7 ->        ; t4*IRF7                                           ;
        E6a: -> IRF7P       ; H*k51*IRF7                                        ;
        E6b: IRF7P ->       ; t5*IRF7P                                          ;
        
        // Conversion factors
        s_t = {s_to_mcs / 60.0 / 60.0} ;
        
        //Parameters
        k11 = {IFNInputs.k11} * s_t ;
        k12 = {IFNInputs.k12} * s_t ;
        k13 = {IFNInputs.k13}       ;
        k14 = {IFNInputs.k14} * s_t ;
        k21 = {IFNInputs.k21} * s_t ;
        k31 = {IFNInputs.k31} * s_t ;
        k32 = {IFNInputs.k32}       ;
        k33 = {IFNInputs.k33}       ;
        t3  = {IFNInputs.t3}  * s_t ;
        k41 = {IFNInputs.k41} * s_t ;
        k42 = {IFNInputs.k42} * s_t ;
        t4  = {IFNInputs.t4}  * s_t ;
        k51 = {IFNInputs.k51} * s_t ;
        t5  = {IFNInputs.t5}  * s_t ;
        n   = {IFNInputs.n}         ;
        RIGI = {IFNInputs.RIGI}     ;

        // Inputs
        H    = {_H}      ;
        IFNe = {_IFNe}   ;
        V    = {_V}      ;
    end'''
    return model_string


def viral_replication_model_string(_v_internalization = 0.0, _IFNe = 0.0):
    """
    Antimony model string generator for viral replication model coupled with the intracellular signaling
    model adopted from "Multiscale Model of RNA Virus Replication and Interferon Responses Reveals Factors
    Controlling Plaque Growth Dynamics" Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego,
    James A. Glazier and Jason E. Shoemaker
    :return: Antimony model string
    """
    model_string = f'''viral_replication_model()
        //Equations
        E7a: H ->           ; H*k61*V                     ;
        E8a: -> V           ; H*k71*V/(1.0+k72*IFNe*7E-5) ;
        E8b: V ->           ; k73*V                       ;
        
        //Parameters
        k61 = {IFNInputs.k61} * s_t ;
        k71 = {IFNInputs.k71} * s_t ;
        k72 = {IFNInputs.k72}       ;
        k73 = {IFNInputs.k73} * s_t ;
        
        //Initial Conditions
        V = {_v_internalization}    ; 
        H = 1.0                     ;
    
        //Inputs
        IFNe = {_IFNe}   ;
    '''
    return model_string

ifn_model_vars = ["IFN", "STATP", "IRF7", "IRF7P"]
viral_replication_model_vars = ["H", "V"]

class IFNDataSteppable(ViralInfectionVTMSteppableBasePy):
    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self.ifn_data_win = None
        self.ifn_data_path = None
        self.ifn_data = dict()

        self.plot_ifn_data = IFNInputs.plot_ifn_data_freq > 0
        self.write_ifn_data = IFNInputs.write_ifn_data_freq > 0

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

    def start(self):
        if self.plot_ifn_data:
            self.ifn_data_win = self.add_new_plot_window(title='IFN Signaling',
                                                          x_axis_title='Time (hrs)',
                                                          y_axis_title='Variables', x_scale_type='linear',
                                                          y_scale_type='log',
                                                          grid=False,
                                                          config_options={'legend': True})

            colors = ['blue', 'red', 'green', 'yellow', 'white', 'cyan', 'purple', 'magenta', 'darkgreen']
            for x in range(len(ifn_model_vars)):
                self.ifn_data_win.add_plot(ifn_model_vars[x], style='Dots', color=colors[x], size=5)

        if self.write_ifn_data:
            from pathlib import Path
            self.ifn_data_path = Path(self.output_dir).joinpath('ifn_data.dat')
            with open(self.ifn_data_path, 'w'):
                pass

    def step(self, mcs):
        if self.simdata_steppable is None:
            self.simdata_steppable = self.shared_steppable_vars[ViralInfectionVTMLib.simdata_steppable_key]

        plot_ifn_data = self.plot_ifn_data and mcs % IFNInputs.plot_ifn_data_freq == 0
        write_ifn_data = self.write_ifn_data and mcs % IFNInputs.write_ifn_data_freq == 0

        #TODO: Modify to plot and output AVG, MAX and MIN
        vrm_tracked_cell = self.simdata_steppable.vrm_tracked_cell
        if vrm_tracked_cell is not None and (plot_ifn_data or write_ifn_data):
            try:
                cell_sbml = getattr(vrm_tracked_cell.sbml, self.vr_model_name)
                if plot_ifn_data:
                    for x in ifn_model_vars:
                        self.ifn_data_win.add_data_point(x, mcs * s_to_mcs / 60.0 / 60.0, cell_sbml[x])

                if write_ifn_data:
                    self.ifn_data[mcs * s_to_mcs / 60.0 / 60.0] = [vrm_tracked_cell.id]
                    for x in ifn_model_vars:
                        self.ifn_data[mcs].append(cell_sbml[x])
            except KeyError:
                pass

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
        if self.write_ifn_data and self.ifn_data:
            with open(self.ifn_data_path, 'a') as fout:
                fout.write(SimDataSteppable.data_output_string(self, self.ifn_data))
                self.ifn_data.clear()