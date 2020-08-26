import sys
import os
from cc3d.core.PySteppables import *

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
from ViralInfectionVTMModelInputs import s_to_mcs, vr_step_size
import ViralInfectionVTMLib
from ViralInfectionVTMSteppables import SimDataSteppable

from .DrugDosingInputs import *

drug_dosing_model_key = "drug_dose_steppable"

days_2_mcs = s_to_mcs / 60 / 60 / 24

class DrugDosingModelSteppable(SteppableBasePy):
    """
    Implements drug dosing regimen
    """

    def __init__(self, frequency=1):
        super().__init__(self, frequency)

        self.drug_dosing_model_key = drug_dosing_model_key

        self.plot_ddm_data = plot_ddm_data_freq > 0
        self.write_ddm_data = write_ddm_data_freq > 0

    def set_drug_model_string(self, _init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4,
                              _k0_rate, _d0_rate, _k1_rate, _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate,
                              _d4_rate, _first_dose, _initial_dose, _dose_interval, _dose):
        """
        Antimony model string generator for this steppable.
        To change parameters do so on the DrugDosingInputs
        :param
        """

        dosingmodel_str = '''

        model dosingmodel()

        // Simple cascade model of bioiavailability with multiple metabolites
        // linear clearance at each stage
        // All times measured in Days

        J0: Drug -> Available1 ; k0*Drug ; //Distribution and bioavailability of drug after dosing
        J0A: Drug -> ; d0*Drug ; // Clearance of drug before bioavailability
        J1: Available1 -> Available2 ; k1*Available1 ; // Metabolism of drug into metabolite 2
        J1A: Available1 -> ; d1*Available1 ; // Clearance of drug after bioavailability
        J2: Available2 -> Available3 ; k2*Available2 ; // Metabolism of drug into metabolite 3
        J2A: Available2 -> ; d2*Available2 ; // Clearance of metabolite 2 
        J3: Available3 -> Available4 ; k3*Available3 ; // Metabolism of drug into metabolite 4
        J3A: Available3 -> ; d3*Available3 ; // Clearance of metabolite 3 
        J4A: Available4 -> ; d4*Available4 ; // Clearance of metabolite 4

        //Initial values
        Drug = {} ; 
        Available1 = {};
        Available2 = {};
        Available3 = {};
        Available4 = {};

        k0 = {}; // bioavailability rate, units /day
        d0 = {} ; // clearance time, units /day 
        k1 = {} ; // metabolism of primary drug rate, units /day
        d1 = {} ; // clearance time, units /day = 4 hours
        k2 = {} ; // metabolism of secondary product, units /day
        d2 = {} ; // clearance time, units /day = 4 hours
        k3 = {} ; // metabolism of tertiary product, units /day
        d3 = {} ; // clearance time, units /day = 4 hours
        d4 = {} ; // clearance time, units /day = 4 hours

        first_dose={} ; // time of first dose in days
        initial_dose = {} ; // initial dose (arbitrary amount)
        dose_interval = {} ; // time interval between doses in days
        dose = {} ; //dose of subsequent treatments

        E1: at (time>first_dose): Drug=Drug+initial_dose ;
        E2: at ( (time-first_dose > 0) && sin((((time-first_dose)/dose_interval))*2*pi)>0): Drug=Drug+dose
        end
        '''.format(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4, _k0_rate, _d0_rate, _k1_rate,
                   _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate, _d4_rate, _first_dose, _initial_dose,
                   _dose_interval, _dose)

        drug_dosig_model_vars = ["Drug", "Available1", "Available2", "Available3", "Available4"]

        return dosingmodel_str, drug_dosig_model_vars

    def start(self):

        # set model string
        self.drug_model_string, self.ddm_vars = self.set_drug_model_string(Drug, Available1, Available2, Available3,
                                                                           Available4,
                                                                           k0, d0, k1, d1, k2, d2, k3, d3, d4,
                                                                           first_dose,
                                                                           initial_dose, dose_interval, dose)

        # init sbml
        self.add_free_floating_antimony(model_string=self.drug_model_string, step_size=days_2_mcs,
                                        model_name='drug_dosing_model')

        # init plots
        if self.plot_ddm_data:
            self.ddm_data_win = self.add_new_plot_window(title='Drug dosing model',
                                                         x_axis_title='Time (hours)',
                                                         y_axis_title='Variables',
                                                         x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True,
                                                         config_options={'legend': True})
            colors = ['blue', 'red', 'green', 'yellow', 'white']
            for c, var in zip(colors, self.ddm_vars):
                self.ddm_data_win.add_plot(var, style='Dots', color=c, size=5)

        # todo: init data save
        # save data

        # Post reference to self
        self.shared_steppable_vars[self.drug_dosing_model_key] = self

    def step(self, mcs):

        if self.plot_ddm_data and mcs % plot_ddm_data_freq == 0:
            [self.ddm_data_win.add_data_point(x, s_to_mcs * mcs / 60 / 60, self.sbml.drug_dosing_model[x])
             for x in self.ddm_vars]

        # todo: modify rmax -> n=2 diminishing hill function
        # todo plot rmax
        # todo do map investigation of max value of avail4 to EC50; ie max(avail4) = [.25, .5, .75, 1, 1.5, 2, 5] EC50

        self.timestep_sbml()

    def on_stop(self):
        self.finish()

    def finish(self):
        # todo: call flush function
        return
