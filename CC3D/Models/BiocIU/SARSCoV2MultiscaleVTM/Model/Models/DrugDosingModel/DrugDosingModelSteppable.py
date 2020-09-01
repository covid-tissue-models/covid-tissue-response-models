import sys
import os
from cc3d.core.PySteppables import *

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
from ViralInfectionVTMModelInputs import s_to_mcs, vr_step_size, replicating_rate, kon, koff, \
    initial_unbound_receptors, hill_coeff_uptake_pr, rate_coeff_uptake_pr, max_ck_secrete_infect, unpacking_rate, \
    r_half, translating_rate, packing_rate, secretion_rate
import ViralInfectionVTMLib
from ViralInfectionVTMSteppableBasePy import *
# ViralInfectionVTMSteppableBasePy.vr_model_name
# from ViralInfectionVTMSteppableBasePy import vr_model_name

from ViralInfectionVTMSteppables import SimDataSteppable
from nCoVToolkit import nCoVUtils

from .DrugDosingInputs import *

drug_dosing_model_key = "drug_dose_steppable"

days_2_mcs = s_to_mcs / 60 / 60 / 24

'''
with the default parameters (k0 = 100.0; d0 = 1.0; k1 = 25.0; d1 = 6.0; k2 = 25.0; d2 = 6.0; k3 = 25.0; d3 = 6.0; 
d4 = 6.0) max(Available4) is a linear function of dose following:

max(Available4) ~= 4.14360796e-01 x dose -1.65564741e-08
'''


class DrugDosingModelSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements drug dosing regimen
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        self.drug_dosing_model_key = drug_dosing_model_key

        self.plot_ddm_data = plot_ddm_data_freq > 0
        self.write_ddm_data = write_ddm_data_freq > 0

        self.max_avail4 = 4.14360796e-01 * dose  # see comment just before steppable definition

        self.vr_model_name = ViralInfectionVTMLib.vr_model_name

        self.__flush_counter = 1

        if self.write_ddm_data:
            self.data_files = {'ddm_data': 'ddm_data.dat', 'ddm_rmax_data': 'ddm_rmax_data.dat'}
            self.ddm_data = {'ddm_data': {}, 'ddm_rmax_data': {}}

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

    def init_plots(self):
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

        self.rmax_data_win = self.add_new_plot_window(title='r_max vs Time',
                                                      x_axis_title='Time (hours)',
                                                      y_axis_title='r_max',
                                                      x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True,
                                                      config_options={'legend': True})

        self.rmax_data_win.add_plot('rmax', style='Dots', color='red', size=5)

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
        self.rmax = self.get_rmax(self.sbml.drug_dosing_model['Available4'])

        # replace viral uptake function
        vim_steppable = self.shared_steppable_vars[ViralInfectionVTMLib.vim_steppable_key]

        vim_steppable.do_cell_internalization = self.do_cell_internalization_w_rmax

        # init plots
        if self.plot_ddm_data:
            self.init_plots()

        # init save data
        if self.write_ddm_data:
            from pathlib import Path
            for key, rel_path in self.data_files.items():
                self.data_files[key] = Path(self.output_dir).joinpath(rel_path)
                with open(self.data_files[key], 'w'):
                    pass

        # Post reference to self
        self.shared_steppable_vars[self.drug_dosing_model_key] = self

    def get_rmax(self, avail4):
        return (1 - nCoVUtils.hill_equation(avail4, rel_avail4_EC50 * self.max_avail4, 2)) * replicating_rate

    def step(self, mcs):
        self.rmax = self.get_rmax(self.sbml.drug_dosing_model['Available4'])
        # print(rmax)
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING):
            vr_model = getattr(cell.sbml, self.vr_model_name)
            vr_model.replicating_rate = self.rmax

        if self.plot_ddm_data and mcs % plot_ddm_data_freq == 0:
            [self.ddm_data_win.add_data_point(x, s_to_mcs * mcs / 60 / 60, self.sbml.drug_dosing_model[x])
             for x in self.ddm_vars]
            if mcs > first_dose / days_2_mcs:
                self.rmax_data_win.add_data_point('rmax', s_to_mcs * mcs / 60 / 60, self.rmax)

        if self.write_ddm_data and mcs % write_ddm_data_freq == 0:
            self.ddm_data['ddm_rmax_data'][mcs] = [self.rmax]

            self.ddm_data['ddm_data'][mcs] = [self.sbml.drug_dosing_model[x] for x in self.ddm_vars]

        # todo do map investigation of max value of avail4 to EC50; ie max(avail4) = [.25, .5, .75, 1, 1.5, 2, 5] EC50

        if mcs >= int(self.simulator.getNumSteps() / 4 * self.__flush_counter):
            self.flush_stored_outputs()
            self.__flush_counter += 1

        self.timestep_sbml()

    def do_cell_internalization_w_rmax(self, cell, viral_amount_com):
        if cell.dict['Receptors'] == 0:
            return False, 0.0

        _k = kon * cell.volume / koff
        diss_coeff_uptake_pr = (initial_unbound_receptors / 2.0 / _k / cell.dict['Receptors']) ** \
                               (1.0 / hill_coeff_uptake_pr)
        uptake_probability = nCoVUtils.hill_equation(viral_amount_com,
                                                     diss_coeff_uptake_pr,
                                                     hill_coeff_uptake_pr)

        cell_does_uptake = np.random.rand() < uptake_probability
        uptake_amount = s_to_mcs / rate_coeff_uptake_pr * uptake_probability

        if cell_does_uptake and cell.type == self.UNINFECTED:
            cell.type = self.INFECTED
            cell.dict['ck_production'] = max_ck_secrete_infect
            self.load_viral_replication_model(cell=cell, vr_step_size=vr_step_size,
                                              unpacking_rate=unpacking_rate,
                                              # replicating_rate=replicating_rate,
                                              replicating_rate=self.rmax,
                                              r_half=r_half,
                                              translating_rate=translating_rate,
                                              packing_rate=packing_rate,
                                              secretion_rate=secretion_rate)

        return cell_does_uptake, uptake_amount

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage
        :return: None
        """
        # Each tuple contains the necessary information for writing a set of data to file
        #   1. Boolean for whether we're writing to file at all
        #   2. The path to write the data to
        #   3. The data to write
        out_info = [(self.write_ddm_data, self.data_files[x], self.ddm_data[x]) for x in self.ddm_data.keys()]

        for write_data, data_path, data in out_info:
            if write_data:
                with open(data_path, 'a') as fout:
                    fout.write(SimDataSteppable.data_output_string(self, data))
                    data.clear()

    def on_stop(self):
        self.finish()

    def finish(self):
        self.flush_stored_outputs()
