# Model effect of Drug Dosing in viral replication
# Written by J. F. Gianlupi, M.S.
# #todo write some more
# Model parameters are specified in DrugDosingInputs.py
#
# DrugDosingModelSteppable
#   Description: implements drug dosing and viral replication rate reduction
#   Usage:
#       In ViralInfectionVTM.py, add the following
#
#           from Models.DrugDosingModel.DrugDosingModelSteppable import DrugDosingModelSteppable
#           CompuCellSetup.register_steppable(steppable=DrugDosingModelSteppable(frequency=1))


import sys
import os
from cc3d.core.PySteppables import *

sys.path.append(os.path.join(os.environ["ViralInfectionVTM"], "Simulation"))
sys.path.append(os.environ["ViralInfectionVTM"])
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

from BatchRun import BatchRunLib

drug_dosing_model_key = "drug_dose_steppable"

days_2_mcs = s_to_mcs / 60 / 60 / 24

# todo: make the target of the drug a model input; i.e. pass the name of the variable to be affected by the drug as
#  input

'''
with the default parameters (k0 = 100.0; d0 = 1.0; k1 = 25.0; d1 = 6.0; k2 = 25.0; d2 = 6.0; k3 = 25.0; d3 = 6.0; 
d4 = 6.0) max(Available4) is a linear function of dose following:

max(Available4) ~= 4.14360796e-01 x dose -1.65564741e-08

with the fitted pk (k0 = 10.0; d0 = 16.635; k1 = 1.0; d1 = 8.317; k2 = 989.6; d2 = 8.317; k3 = 158.4; d3 = 0.693; 
d4 = 0.693) max(Available4) is a linear function of dose following:

max(Available4) ~= 2.32417475e-01 x dose + 1.59151098e-08

'''


# @staticmethod
def set_default_ddm_string(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4,
                           _k0_rate, _d0_rate, _k1_rate, _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate,
                           _d4_rate, _first_dose, _initial_dose, _dose_interval, _dose, _eot):
    """
    Antimony model string generator for this steppable.
    To change parameters do so on the DrugDosingInputs. Parameter descriptions are also in DrugDosingInputs
    :param
    """

    dosingmodel_str = '''
    model dosingmodel()
    //Time is in days!
    
    //infusion
    
    -> Dpls; switch * infusion_amount / one_our // switch = (0,1) to turn on or off, infusion happens over 1h
    
    //flow from plasma 
    
    Dpls -> ; kE0 * Dpls // elimination
    
    Dpls -> Dperi ; kp * Dpls // to periphery 
    
    Dpls -> Dlung ; k0 * Dpls
    
    // flow from periphery
    
    Dperi -> Dpls ; kpp * Dpls // to plasma 
    
    // Drug reactions / flow in lung
    
    Dlung -> Dpls ; k0 * Dlung
    
    Dlung -> Mala ; k12 * Dlung
    
    Dlung -> ; kE1 * Dlung
    
    // Mala reactions
    
    Mala -> Mnmp ; k23 * Mala
    
    Mala -> ; kE2 * Mala
    
    //Mnmp reactions
    
    Mnmp -> Mntp ; k34 * Mnmp
    Mnmp ->  ; kE3 * Mnmp
    
    // Mntp reaction
    
    Mntp -> ; kE4 * Mntp
    
    //parameters
    // initial conditions 
    
    Dpls = {}
    
    Dperi = {}
    
    Dlung = {}
    
    Mala = {}
    
    Mnmp = {}
    
    Mntp = {}
    
    
    //utils
    switch = 0 //turns infusion on/off
    curr_infu_start = 0 // tracks when current infusion started
    
    // rates
    
    kp = 0.41195
    
    kpp = 0.36502

    k0 = 6.3335
    
    k12 = 1.2248
    
    k23 = 372.61
    
    k34 = 181.64

    kE0 = 20.253
    
    kE1 = 7.81
    
    kE2 = 6.0801
    
    kE3 = 0.97259
    
    kE4 = 0.83115

    //constants
    infusion_amount = 1
    
    dose_interval = 1 // time interval between doses in days
    
    dose_end = 9999 // end of treatment day
    
    one_our = 1/24 
    
    first_dose = 0 // time of first dose in days
    
    // events
    
    E1: at (time - first_dose > 0): switch = 1, curr_infu_start = time ; // starts the first infusion
    E2: at ( (time-first_dose > dose_interval) && (time < dose_end) && sin((((time-first_dose)/dose_interval))*2*pi)>0): switch = 1, curr_infu_start = time; // starts the subsequent infusions
    E3: at (time - (one_our + curr_infu_start) > 0): switch = 0 ; // turns infusion off
    
    end
'''.format(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4, _k0_rate, _d0_rate, _k1_rate,
                       _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate, _d4_rate, _first_dose, _initial_dose,
                       _dose_interval, _dose, _eot)

    drug_dosig_model_vars = ["Drug", "Available1", "Available2", "Available3", "Available4"]

    return dosingmodel_str, drug_dosig_model_vars


def full_ddm_for_testing(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4,
                         _k0_rate, _d0_rate, _k1_rate, _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate,
                         _d4_rate, _first_dose, _initial_dose, _dose_interval, _dose, _eot):
    """
    Antimony model string generator for this steppable.
    To change parameters do so on the DrugDosingInputs. Parameter descriptions are also in DrugDosingInputs
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
            dose_end = {} // end of treatment day

            E1: at (time - first_dose > 0): Drug=Drug+initial_dose ;
            E2: at ( (time-first_dose > dose_interval) && (time < dose_end) && sin((((time-first_dose)/dose_interval))*2*pi)>0): Drug=Drug+dose
            end
            '''.format(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4, _k0_rate, _d0_rate, _k1_rate,
                       _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate, _d4_rate, _first_dose, _initial_dose,
                       _dose_interval, _dose, _eot)

    drug_dosig_model_vars = ["Drug", "Available1", "Available2", "Available3", "Available4"]

    return dosingmodel_str, drug_dosig_model_vars


def set_cst_drug_ddm_string(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4,
                            _k0_rate, _d0_rate, _k1_rate, _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate,
                            _d4_rate, _first_dose, _initial_dose, _dose_interval, _dose, _eot):
    """
    Antimony model string generator for this steppable.
    To change parameters do so on the DrugDosingInputs. Parameter descriptions are also in DrugDosingInputs
    :param
    """

    dosingmodel_str = '''

            model dosingmodel()

            // Simple cascade model of bioiavailability with multiple metabolites
            // linear clearance at each stage
            // All times measured in Days

           // J0: Drug -> Available1 ; k0*Drug ; //Distribution and bioavailability of drug after dosing
            J0A: Drug -> ; d0*Drug ; // Clearance of drug before bioavailability
            J1: Available1 -> Available2 ; k1*Available1 ; // Metabolism of drug into metabolite 2
            J1A: Available1 -> ; d1*Available1 ; // Clearance of drug after bioavailability
            J2: Available2 -> Available3 ; k2*Available2 ; // Metabolism of drug into metabolite 3
            J2A: Available2 -> ; d2*Available2 ; // Clearance of metabolite 2 
            J3: Available3 -> Available4 ; k3*Available3 ; // Metabolism of drug into metabolite 4
            J3A: Available3 -> ; d3*Available3 ; // Clearance of metabolite 3 
            J4A: Available4 -> ; d4*Available4 ; // Clearance of metabolite 4

            //Initial values
            dummy = {} ; 
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
            dose_end = {} // end of treatment day

            const Drug := initial_dose;

            //E1: at (time - first_dose > 0): Drug=Drug+initial_dose ;
            //E2: at ( (time-first_dose > dose_interval) && (time < dose_end) && sin((((time-first_dose)/dose_interval))*2*pi)>0): Drug=Drug+dose
            end
            '''.format(_init_drug, _init_avail1, _init_avail2, _init_avail3, _init_avail4, _k0_rate, _d0_rate, _k1_rate,
                       _d1_rate, _k2_rate, _d2_rate, _k3_rate, _d3_rate, _d4_rate, _first_dose, _initial_dose,
                       _dose_interval, _dose, _eot)

    drug_dosig_model_vars = ["Drug", "Available1", "Available2", "Available3", "Available4"]

    return dosingmodel_str, drug_dosig_model_vars


def set_cell_drug_metabolization(_init_avail1, _init_avail2, _init_avail3, _init_avail4, _k0_rate, _k1_rate, _d1_rate,
                                 _k2_rate, _d2_rate, _k3_rate, _d3_rate, _d4_rate):
    """

    Antimony model string generator for drug metabolization in cells.
    To change parameters do so on the DrugDosingInputs. Parameter descriptions are also in DrugDosingInputs


    :param _k0_rate:
    :param _init_avail1:
    :param _init_avail2:
    :param _init_avail3:
    :param _init_avail4:
    :param _k1_rate:
    :param _d1_rate:
    :param _k2_rate:
    :param _d2_rate:
    :param _k3_rate:
    :param _d3_rate:
    :param _d4_rate:
    :return:
    """

    metabolization_str = '''

                model metabolizationmodel()

                // Simple cascade model of bioiavailability with multiple metabolites
                // linear clearance at each stage
                // All times measured in Days

                
                J1: Available1 -> Available2 ; k1*Available1 ; // Metabolism of drug into metabolite 2
                J1A: Available1 -> ; d1*Available1 ; // Clearance of drug after bioavailability
                J2: Available2 -> Available3 ; k2*Available2 ; // Metabolism of drug into metabolite 3
                J2A: Available2 -> ; d2*Available2 ; // Clearance of metabolite 2 
                J3: Available3 -> Available4 ; k3*Available3 ; // Metabolism of drug into metabolite 4
                J3A: Available3 -> ; d3*Available3 ; // Clearance of metabolite 3 
                J4A: Available4 -> ; d4*Available4 ; // Clearance of metabolite 4

                //Initial values
                 
                Available1 = {};
                Available2 = {};
                Available3 = {};
                Available4 = {};
                
                k0 = {}; // bioavailability rate, units /day
                k1 = {} ; // metabolism of primary drug rate, units /day
                d1 = {} ; // clearance time, units /day = 4 hours
                k2 = {} ; // metabolism of secondary product, units /day
                d2 = {} ; // clearance time, units /day = 4 hours
                k3 = {} ; // metabolism of tertiary product, units /day
                d3 = {} ; // clearance time, units /day = 4 hours
                d4 = {} ; // clearance time, units /day = 4 hours

                end
                '''.format(_init_avail1, _init_avail2, _init_avail3, _init_avail4, _k0_rate, _k1_rate, _d1_rate,
                           _k2_rate, _d2_rate, _k3_rate, _d3_rate, _d4_rate)

    return metabolization_str


class DrugDosingModelSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements drug dosing regimen
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        import Models.DrugDosingModel.DrugDosingInputs as DrugDosingInputs
        BatchRunLib.apply_external_multipliers(__name__, DrugDosingInputs)
        self.drug_dosing_model_key = drug_dosing_model_key

        if constant_drug_concentration:
            self.set_drug_model_string = set_cst_drug_ddm_string
        else:
            self.set_drug_model_string = set_default_ddm_string

        self.set_control_model_string = full_ddm_for_testing

        self.plot_ddm_data = plot_ddm_data_freq > 0
        self.write_ddm_data = write_ddm_data_freq > 0

        self.max_avail4 = 2.32417475e-01 * dose  # see comment just before steppable definition

        if auto_ec50:
            self.hill_k = self.max_avail4 * rel_avail4_EC50
        else:
            self.hill_k = ec50

        self.drug_model_string = None

        self.ddm_vars = None

        self.drug_metabolization_string = None

        self.control_string = None

        self.rmax = None

        self.vr_model_name = ViralInfectionVTMLib.vr_model_name

        self.ddm_rr = None

        self.control_rr = None

    @staticmethod
    def get_roadrunner_for_single_antimony(model):
        """
        :type model: str name of the model
        :param model:
        :return:
        """
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == model:
                return rr
        return None

    @staticmethod
    def get_sbml_simulator_for_cell(model_name: str, cell: object = None) -> Union[object, None]:
        """
        Returns a reference to RoadRunnerPy or None
        :param model_name: model name
        :param cell: CellG cell object
        :return {instance of RoadRunnerPy} or {None}:
        """
        try:
            dict_attrib = CompuCell.getPyAttrib(cell)
            return dict_attrib['SBMLSolver'][model_name]
        except LookupError:
            return None

    def timestep_cell_sbml(self, model_name: str, cell: object = None):
        if not cell:
            return
        rr = self.get_sbml_simulator_for_cell(model_name, cell)
        rr.timestep()

    def start(self):

        # set model string
        self.drug_model_string, self.ddm_vars = self.set_drug_model_string(Drug, Available1, Available2, Available3,
                                                                           Available4,
                                                                           k0, d0, k1, d1, k2, d2, k3, d3, d4,
                                                                           first_dose,
                                                                           initial_dose, dose_interval, dose, dose_end)
        self.control_string, _ = self.set_control_model_string(Drug, Available1, Available2, Available3,
                                                               Available4,
                                                               k0, d0, k1, d1, k2, d2, k3, d3, d4,
                                                               first_dose,
                                                               initial_dose, dose_interval, dose, dose_end)

        # init sbml
        self.add_free_floating_antimony(model_string=self.drug_model_string, step_size=days_2_mcs,
                                        model_name='drug_dosing_model')
        self.ddm_rr = self.get_roadrunner_for_single_antimony('drug_dosing_model')

        self.add_free_floating_antimony(model_string=self.control_string, step_size=days_2_mcs,
                                        model_name='drug_dosing_control')
        self.control_rr = self.get_roadrunner_for_single_antimony('drug_dosing_control')

        self.drug_metabolization_string = set_cell_drug_metabolization(Available1, Available2, Available3, Available4,
                                                                       k0, k1, d1, k2, d2, k3, d3, d4)

        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
            self.add_antimony_to_cell(model_string=self.drug_metabolization_string, model_name='drug_metabolization',
                                      cell=cell, step_size=days_2_mcs)

        if prophylactic_treatment:
            # to be able to write the data from prophylaxis I put the prophylactic code in the
            # data steppable. May not be elegant but it works
            # this DOES MEAN that if the write step is not included prophylaxis won't work
            pass

        if sanity_run:
            self.rmax = replicating_rate
        else:
            self.rmax = self.get_rmax(self.sbml.drug_dosing_model['Available4'])

        self.shared_steppable_vars['rmax'] = self.rmax

        # replace viral uptake function
        vim_steppable = self.shared_steppable_vars[ViralInfectionVTMLib.vim_steppable_key]

        vim_steppable.do_cell_internalization = self.do_cell_internalization_changing_rmax

        # Post reference to self
        self.shared_steppable_vars[self.drug_dosing_model_key] = self

    def get_rmax(self, avail4):
        return (1 - nCoVUtils.hill_equation(avail4, self.hill_k, 2)) * replicating_rate

    def get_all_uptakes_scalar_prodrug(self):
        """

        :return:
        """

        total = 0
        uptakes = []
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
            rate = days_2_mcs * cell.sbml.drug_metabolization['k0']  # k0 is in units of /day!!!!!!!!!!!!!
            u = rate * self.sbml.drug_dosing_model['Drug'] / len(self.cell_list_by_type(
                self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED))
            total += u
            uptakes.append((cell.id, u))

        return uptakes, total

    def do_prodrug_metabolization(self):
        """

        :return:
        """
        if self.sbml.drug_dosing_model['Drug'] <= 0:
            print(self.sbml.drug_dosing_model['Drug'])
            if self.sbml.drug_dosing_model['Drug'] < 0:
                self.sbml.drug_dosing_model['Drug'] = 0
            return 0
        if not diffusing_drug:

            # for non-diffusing drug (aka, an scalar) the calculation of uptake goes as follows:
            # uptake_total = sum(uptake_cell) = sum( rate_cell * drug_cell); as non-diffusing drug_cell is the same for
            # all cells, call it pc
            # uptake_total = pc * sum(rate_cell); say now that all cells metabolize at the same rate, kc
            # uptake_total = pc * kc * Nc; Nc being the number of cells
            # It follows, then
            # uptake_total = Nc * uptake_cell.
            # We know the total uptake from the ODE,
            # uptake_total = k_ode * drug_total, so
            # Nc * uptake_cell = uptake_total = k_ode * drug_total
            # kc * pc = uptake_cell = k_ode * drug_total / Nc

            uptakes, total = self.get_all_uptakes_scalar_prodrug()

            if total < self.sbml.drug_dosing_model['Drug']:
                for cid, u in uptakes:
                    cell = self.fetch_cell_by_id(cid)
                    cell.sbml.drug_metabolization['Available1'] += u
                    # print(cell.id, cell.sbml.drug_metabolization['Available1'])
                    self.sbml.drug_dosing_model['Drug'] -= u

            else:
                total = self.sbml.drug_dosing_model['Drug']
                u = self.sbml.drug_dosing_model['Drug'] / len(self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING,
                                                                                     self.UNINFECTED))
                # print('equal uptake = ', u)
                for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                    cell.sbml.drug_metabolization['Available1'] += u
                    self.sbml.drug_dosing_model['Drug'] -= u

        else:
            # regular cc3d uptake
            secreter = self.get_field_secretor("prodrug")
            total = 0
            for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                rate = days_2_mcs * cell.sbml.drug_metabolization['k0']
                uptake = secreter.uptakeInsideCellTotalCount(cell, 9e99, rate)  # uptake at rate 'rate' without max
                # value for uptake (9e99)
                cell.sbml.drug_metabolization['Available1'] += abs(uptake.tot_amount)
                total += abs(uptake.tot_amount)
        print('result:', self.sbml.drug_dosing_model['Drug'], total)
        print('control:', self.sbml.drug_dosing_control['Drug'], days_2_mcs * self.sbml.drug_dosing_control['J0'])
        print('result/control:', self.sbml.drug_dosing_model['Drug'] / self.sbml.drug_dosing_control['Drug'],
              total / (days_2_mcs * self.sbml.drug_dosing_control['J0']))
        return total

    def step(self, mcs):

        # the rate k0 is to go
        # from plasma to the lung and back. So how much enters the system is k0 (Dplasma - Dlung).
        # Then the cells convert it.
        # What we want then is to have the chemical be deposited in the epithelial cells. We also want it to not
        # diffuse to medium (if possible). Cells read the full value in their region and uptake k12*Dlung(uptaken)
        # 0. Get the amount leaving back to plasma
        # 1. get amount of prodrug entering the system, k0*D, add it (D = Dplasma - Dlung).
        # 1.1. deposit k0*D uniformely as diffusing concentration
        # 1.2. the amount leaving the system is modeled by the decay. Decay will be both k0 and kE1, \gamma = k0 + kE1;
        # so need to calculate k0*Dlung to add it back to the sbml
        # 1.note I'll have the prodrug as a global for the beginning of implementation, change later
        #
        # 2. The epithelial cells uptake in their whole domain.
        # 2.1. They uptake k12*Dlung(over cell) and that is added to their sbml
        #
        # CELLULARIZATION NOTE!!!
        # For drug uptake need to properly distribute it. Say the total dose given is 1, I need to divide the total dose
        # by the person's weight, then multiply it by the weight of the patch, then divide by the number of cells the
        # patch initially had.
        # !!!!!!!!!!!!!!!!!!!!NOTE ON THE NOTE!!!!!!!!!!!!!!!!!!!!
        # It's wrong, overthinking on my part (we are using concentrations, so all good)
        #
        #
        #
        #
        # diffusion of remdesivir: very fast, mol weigh of 602.585. See bose-einstein for upper limit. treat it as a
        # small molecule.

        self.ddm_rr.timestep()
        self.control_rr.timestep()
        remdesivir_upt_tot = self.do_prodrug_metabolization()
        self.shared_steppable_vars['remdesivir_upt_tot'] = remdesivir_upt_tot
        if not sanity_run:
            # self.rmax = self.get_rmax(self.sbml.drug_dosing_model['Available4'])
            # self.shared_steppable_vars['rmax'] = self.rmax

            for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                self.timestep_cell_sbml('drug_metabolization', cell)

                cell.dict['rmax'] = self.get_rmax(cell.sbml.drug_metabolization['Available4'])
                if cell.type != self.UNINFECTED:
                    vr_model = getattr(cell.sbml, self.vr_model_name)
                    vr_model.replicating_rate = cell.dict['rmax']
                    # print(vr_model.replicating_rate)

                time = cell.sbml.drug_metabolization['Time']
            print('time', time, self.sbml.drug_dosing_model['Time'])

        # for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING):
        #     vr_model = getattr(cell.sbml, self.vr_model_name)
        #     vr_model.replicating_rate = self.rmax

            # ViralInfectionVTMLib.step_sbml_model_cell(cell=cell)

    def get_rna_array(self):
        return np.array([cell.dict['Replicating'] for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING,
                                                                                     self.UNINFECTED, self.DYING)])

    def do_cell_internalization_changing_rmax(self, cell, viral_amount_com):
        # WARNING!! OVERWRITES FUNCTION OF MAIN MODEL
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

    def finish(self):
        pass


class DrugDosingDataFieldsPlots(ViralInfectionVTMSteppableBasePy):
    """
    Responsible for plots, extra fields and data handling
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        # import Models.DrugDosingModel.DrugDosingInputs as DrugDosingInputs
        self.track_cell_level_scalar_attribute(field_name='internal_viral_RNA', attribute_name='Replicating')

        self.mvars = None

        self.get_rna_array = None

        self.plot_ddm_data = plot_ddm_data_freq > 0
        self.write_ddm_data = write_ddm_data_freq > 0

        self.ddm_data_win = None

        self.ddm_control_plot = None

        self.rmax_data_win = None

        self.total_rna_plot = None
        self.mean_rna_plot = None

        self.__flush_counter = 1

        if self.write_ddm_data:
            self.data_files = {'ddm_data': 'ddm_data.dat', 'ddm_rmax_data': 'ddm_rmax_data.dat',
                               'ddm_tot_RNA_data': 'ddm_tot_RNA_data.dat', 'ddm_mean_RNA_data': 'ddm_mean_RNA_data.dat'}
            self.ddm_data = {'ddm_data': {}, 'ddm_rmax_data': {}, 'ddm_tot_RNA_data': {}, 'ddm_mean_RNA_data': {}}

    def init_plots(self):
        self.ddm_data_win = self.add_new_plot_window(title='Drug dosing model',
                                                     x_axis_title='Time (hours)',
                                                     y_axis_title='Variables',
                                                     x_scale_type='linear',
                                                     y_scale_type='linear',
                                                     grid=True,
                                                     config_options={'legend': True})

        self.ddm_control_plot = self.add_new_plot_window(title='Drug dosing control plot',
                                                         x_axis_title='Time (hours)',
                                                         y_axis_title='Variables',
                                                         x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True,
                                                         config_options={'legend': True})

        colors = ['blue', 'red', 'green', 'yellow', 'white']
        ddm_vars = self.mvars.ddm_vars
        for c, var in zip(colors, ddm_vars):
            self.ddm_data_win.add_plot(var, style='Dots', color=c, size=5)
            self.ddm_control_plot.add_plot(var, style='Dots', color=c, size=5)

        self.rmax_data_win = self.add_new_plot_window(title='Mean r_max vs Time',
                                                      x_axis_title='Time (hours)',
                                                      y_axis_title='r_max',
                                                      x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True,
                                                      config_options={'legend': True})

        self.rmax_data_win.add_plot('rmax', style='Dots', color='red', size=5)
        self.total_rna_plot = self.add_new_plot_window(title='Total internal viral RNA',
                                                       x_axis_title='Time (hours)',
                                                       y_axis_title='Variables',
                                                       x_scale_type='linear',
                                                       y_scale_type='linear',
                                                       grid=True,
                                                       config_options={'legend': True})
        self.total_rna_plot.add_plot('RNA_tot', style='Dots', color='red', size=5)
        self.mean_rna_plot = self.add_new_plot_window(title='Mean internal viral RNA',
                                                      x_axis_title='Time (hours)',
                                                      y_axis_title='Variables',
                                                      x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True,
                                                      config_options={'legend': True})
        self.mean_rna_plot.add_plot('RNA_mean', style='Dots', color='red', size=5)

    def init_writes(self):
        # init save data
        # if self.write_ddm_data:
        from pathlib import Path
        for key, rel_path in self.data_files.items():
            self.data_files[key] = Path(self.output_dir).joinpath(rel_path)
            with open(self.data_files[key], 'w'):
                pass

    def start(self):
        self.mvars = self.shared_steppable_vars[drug_dosing_model_key]  # main ddm class vars
        self.get_rna_array = self.mvars.get_rna_array
        if self.plot_ddm_data:
            self.init_plots()

        if self.write_ddm_data:
            self.init_writes()
        if prophylactic_treatment:
            # from cc3d.CompuCellSetup import persistent_globals as pg
            # for model_name, rr in pg.free_floating_sbml_simulators.items():
            #     if model_name == 'drug_dosing_model':
            #         ddm_rr = rr
            #         break
            number_of_prophylactic_steps = int(prophylactic_time / days_2_mcs)

            ddm_rr = self.shared_steppable_vars[drug_dosing_model_key].ddm_rr
            get_rmax = getattr(DrugDosingModelSteppable, 'get_rmax')
            for i in range(number_of_prophylactic_steps):  # let it run for prophylactic_time days
                # print('time stepping', i)
                ddm_rr.timestep()
                self.shared_steppable_vars['rmax'] = get_rmax(self.mvars, self.sbml.drug_dosing_model['Available4'])
                self.shared_steppable_vars['pre_sim_time'] = number_of_prophylactic_steps - i
                if self.write_ddm_data:
                    self.do_writes(0)
            if self.write_ddm_data:
                self.flush_stored_outputs()
                self.__flush_counter -= 1

    def get_metabolites_in_cell(self, cell):
        # print(cell.sbml.drug_metabolization['Available1'])
        metabolites = [cell.sbml.drug_metabolization[x] for x in self.mvars.ddm_vars[1:]]
        # print(cell.id, l)
        return metabolites

    def get_total_metabolites_in_cells(self):

        m = [[] for x in self.mvars.ddm_vars[1:]]
        # print(m)
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
            cm = self.get_metabolites_in_cell(cell)
            for i in range(len(cm)):
                m[i].append(cm[i])
                # print(m[i])

        return m

    def get_mean_std_rmax(self):
        rmax_list = [cell.dict['rmax'] for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING,
                                                                          self.UNINFECTED)]
        return np.mean(rmax_list), np.std(rmax_list)

    def do_plots(self, mcs):
        """
        :parameter mcs
        :return None
        """

        [self.ddm_control_plot.add_data_point(x, s_to_mcs * mcs / 60 / 60, self.sbml.drug_dosing_control[x])
         for x in self.mvars.ddm_vars]

        self.ddm_data_win.add_data_point('Drug', s_to_mcs * mcs / 60 / 60, self.sbml.drug_dosing_model['Drug'])

        total_mets = self.get_total_metabolites_in_cells()

        for i, x in enumerate(self.mvars.ddm_vars[1:]):
            # print(x, np.sum(total_mets[i]), total_mets[i])
            y = np.sum(total_mets[i])
            # print(y)
            self.ddm_data_win.add_data_point(x, s_to_mcs * mcs / 60 / 60, y)

        if mcs > first_dose / days_2_mcs or constant_drug_concentration:
            mean, _ = self.get_mean_std_rmax()
            self.rmax_data_win.add_data_point('rmax', s_to_mcs * mcs / 60 / 60, mean)

        rna_list = self.get_rna_array()

        self.total_rna_plot.add_data_point('RNA_tot', s_to_mcs * mcs / 60 / 60, np.sum(rna_list))
        self.mean_rna_plot.add_data_point('RNA_mean', s_to_mcs * mcs / 60 / 60, np.mean(rna_list))

    def do_writes(self, mcs):
        if prophylactic_treatment and mcs == 0:
            time = mcs - self.shared_steppable_vars['pre_sim_time']
        else:
            time = mcs
        mean, _ = self.get_mean_std_rmax()
        # self.ddm_data['ddm_rmax_data'][time] = [self.shared_steppable_vars['rmax']]
        self.ddm_data['ddm_rmax_data'][time] = [mean]

        self.ddm_data['ddm_data'][time] = [self.sbml.drug_dosing_model[x] for x in self.mvars.ddm_vars]

        if time >= 0:
            rna_list = self.get_rna_array()

            self.ddm_data['ddm_tot_RNA_data'][mcs] = [np.sum(rna_list)]
            self.ddm_data['ddm_mean_RNA_data'][mcs] = [np.mean(rna_list)]

        if mcs >= int(self.simulator.getNumSteps() / 4 * self.__flush_counter):
            self.flush_stored_outputs()

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
        self.__flush_counter += 1

    def step(self, mcs):

        if self.plot_ddm_data and mcs % plot_ddm_data_freq == 0:
            self.do_plots(mcs)
        if self.write_ddm_data and mcs % write_ddm_data_freq == 0:
            self.do_writes(mcs)

    def on_stop(self):
        self.finish()

    def finish(self):
        if self.write_ddm_data:
            self.flush_stored_outputs()
