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
hour_2_mcs = s_to_mcs / 60 / 60

# @staticmethod
def set_default_ddm_string(_init_drug_plasma, _init_drug_periphery, _init_drug_lung, _init_met_alanine, _init_met_NMP,
                           _init_met_NTP, _double_first_dose, _k_p_rate, _k_p_prime_rate, _k01_rate, _k10_rate,
                           _k12_rate, _k23_rate, _k34_rate, _kE0_rate, _kE1_rate, _kE2_rate, _kE3_rate, _kE4_rate,
                           _infusion_amount, _dose_interval, _eot, _treatment_start):
    """
    Antimony model string generator for this steppable.
    To change parameters do so on the DrugDosingInputs. Parameter descriptions are also in DrugDosingInputs
    :param _treatment_start:
    :param _init_drug_plasma:
    :param _init_drug_periphery:
    :param _init_drug_lung:
    :param _init_met_alanine:
    :param _init_met_NMP:
    :param _init_met_NTP:
    :param _double_first_dose:
    :param _k_p_rate:
    :param _k_p_prime_rate:
    :param _k01_rate:
    :param _k10_rate:
    :param _k12_rate:
    :param _k23_rate:
    :param _k34_rate:
    :param _kE0_rate:
    :param _kE1_rate:
    :param _kE2_rate:
    :param _kE3_rate:
    :param _kE4_rate:
    :param _infusion_amount:
    :param _dose_interval:
    :param _eot:
    :param
    """

    dosingmodel_str = '''
    model dosingmodel()
    //Time is in days!
    
    //infusion    
    J0: -> Dpls; switch * infusion_amount / one_our // switch = (0,1) to turn on or off, infusion happens over 1h    
    //flow from plasma     
    J1: Dpls -> ; kE0 * Dpls // elimination    
    J2: Dpls -> Dperi ; kp * Dpls // to periphery     
    J3: Dpls -> Dlung ; k01 * Dpls
        
    // flow from periphery    
    J4: Dperi -> Dpls ; kpp * Dperi // to plasma     
    
    // Drug reactions / flow in lung    
    J5: Dlung -> Dpls ; k10 * Dlung    
    // J6: Dlung -> Mala ; k12 * Dlung  // this happens in cells  
    J7: Dlung -> ; kE1 * Dlung
    //parameters
    // initial conditions     
    Dpls = {}    
    Dperi = {}    
    Dlung = {}    
 
    //utils
    switch = 0 //turns infusion on/off
    curr_infu_start = 0 // tracks when current infusion started
    double_first_dose = {}
    
    // rates    
    kp = {}    
    kpp = {}
    k01 = {}   
    k10 = {}   
    
    kE0 = {}    
    kE1 = {}    
        
    //constants
    infusion_amount = {}    
    dose_interval = {} // time interval between doses in days    
    dose_end = {} // end of treatment day    
    one_our = 1/24     
    first_dose = {} // time of first dose in days
    
    // events    
    E1: at (time - first_dose > 0): switch = 1*double_first_dose, curr_infu_start = time ; // starts the first infusion
    E2: at ( (time-first_dose > dose_interval) && (time < dose_end) && sin((((time-first_dose)/dose_interval))*2*pi)>0): switch = 1, curr_infu_start = time; // starts the subsequent infusions
    E3: at (time - (one_our + curr_infu_start) > 0): switch = 0 ; // turns infusion off
    
    end
'''.format(_init_drug_plasma, _init_drug_periphery, _init_drug_lung,
           _double_first_dose, _k_p_rate, _k_p_prime_rate, _k01_rate, _k10_rate, _kE0_rate,
           _kE1_rate, _infusion_amount, _dose_interval, _eot, _treatment_start)

    drug_dosig_model_vars = ["Dpls", "Dperi", "Dlung", "Mala", "Mnmp", "Mntp"]

    return dosingmodel_str, drug_dosig_model_vars


def full_ddm_for_testing(_init_drug_plasma, _init_drug_periphery, _init_drug_lung, _init_met_alanine, _init_met_NMP,
                         _init_met_NTP, _double_first_dose, _k_p_rate, _k_p_prime_rate, _k01_rate, _k10_rate, _k12_rate,
                         _k23_rate, _k34_rate, _kE0_rate, _kE1_rate, _kE2_rate, _kE3_rate, _kE4_rate,
                         _infusion_amount, _dose_interval, _eot, _treatment_start):
    dosingmodel_str = '''
    model dosingmodel()
    //Time is in days!
    
    //infusion    
    J0: -> Dpls; switch * infusion_amount / one_our // switch = (0,1) to turn on or off, infusion happens over 1h    
    //flow from plasma     
    J1: Dpls -> ; kE0 * Dpls // elimination    
    J2: Dpls -> Dperi ; kp * Dpls // to periphery     
    J3: Dpls -> Dlung ; k01 * Dpls
        
    // flow from periphery    
    J4: Dperi -> Dpls ; kpp * Dperi // to plasma     
    
    // Drug reactions / flow in lung    
    J5: Dlung -> Dpls ; k10 * Dlung    
    J6: Dlung -> Mala ; k12 * Dlung    
    J7: Dlung -> ; kE1 * Dlung
    
    // Mala reactions    
    J8: Mala -> Mnmp ; k23 * Mala    
    J9: Mala -> ; kE2 * Mala
    
    //Mnmp reactions    
    J10: Mnmp -> Mntp ; k34 * Mnmp
    J11: Mnmp ->  ; kE3 * Mnmp
    
    // Mntp reaction    
    J12: Mntp -> ; kE4 * Mntp
    
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
    double_first_dose = {}
    
    // rates    
    kp = {}    
    kpp = {}
    k01 = {}    
    k10 = {}    
    k12 = {}    
    k23 = {}    
    k34 = {}
    kE0 = {}    
    kE1 = {}    
    kE2 = {}    
    kE3 = {}    
    kE4 = {}
    
    //constants
    infusion_amount = {}    
    dose_interval = {} // time interval between doses in days    
    dose_end = {} // end of treatment day    
    one_our = 1/24     
    first_dose = {} // time of first dose in days
    
    // events    
    E1: at (time - first_dose > 0): switch = 1*double_first_dose, curr_infu_start = time ; // starts the first infusion
    E2: at ( (time-first_dose > dose_interval) && (time < dose_end) && sin((((time-first_dose)/dose_interval))*2*pi)>0): switch = 1, curr_infu_start = time; // starts the subsequent infusions
    E3: at (time - (one_our + curr_infu_start) > 0): switch = 0 ; // turns infusion off
    
    end
'''.format(_init_drug_plasma, _init_drug_periphery, _init_drug_lung, _init_met_alanine, _init_met_NMP, _init_met_NTP,
           _double_first_dose, _k_p_rate, _k_p_prime_rate, _k01_rate, _k10_rate, _k12_rate, _k23_rate, _k34_rate,
           _kE0_rate,
           _kE1_rate, _kE2_rate, _kE3_rate, _kE4_rate, _infusion_amount, _dose_interval, _eot, _treatment_start)

    drug_dosig_model_vars = ["Dpls", "Dperi", "Dlung", "Mala", "Mnmp", "Mntp"]

    return dosingmodel_str, drug_dosig_model_vars


def set_cst_drug_ddm_string(_init_drug_plasma, _init_drug_periphery, _init_drug_lung, _init_met_alanine, _init_met_NMP,
                            _init_met_NTP, _double_first_dose, _k_p_rate, _k_p_prime_rate, _k01_rate, _k10_rate,
                            _k12_rate,
                            _k23_rate, _k34_rate, _kE0_rate, _kE1_rate, _kE2_rate, _kE3_rate, _kE4_rate,
                            _infusion_amount, _dose_interval, _eot, _treatment_start):
    dosingmodel_str = '''
    model dosingmodel()
    //Time is in days!
    
    //infusion    
    //J0: -> Dpls; switch * infusion_amount / one_our // switch = (0,1) to turn on or off, infusion happens over 1h    
    //flow from plasma     
    J1: Dpls -> ; kE0 * Dpls // elimination    
    J2: Dpls -> Dperi ; kp * Dpls // to periphery     
    J3: Dpls -> Dlung ; k01 * Dpls
        
    // flow from periphery    
    J4: Dperi -> Dpls ; kpp * Dperi // to plasma     
    
    // Drug reactions / flow in lung    
    J5: Dlung -> Dpls ; k10 * Dlung    
    J6: Dlung -> Mala ; k12 * Dlung    
    J7: Dlung -> ; kE1 * Dlung
    
    // Mala reactions    
    J8: Mala -> Mnmp ; k23 * Mala    
    J9: Mala -> ; kE2 * Mala
    
    //Mnmp reactions    
    J10: Mnmp -> Mntp ; k34 * Mnmp
    J11: Mnmp ->  ; kE3 * Mnmp
    
    // Mntp reaction    
    J12: Mntp -> ; kE4 * Mntp
    
    //parameters
    // initial conditions     
    Dperi = {}    
    Dlung = {}    
    Mala = {}    
    Mnmp = {}    
    Mntp = {}    
    
    //utils
    switch = 0 //turns infusion on/off
    curr_infu_start = 0 // tracks when current infusion started
    double_first_dose = {}
    
    // rates    
    kp = {}    
    kpp = {}
    k01 = {}    
    k10 = {}    
    k12 = {}    
    k23 = {}    
    k34 = {}
    kE0 = {}    
    kE1 = {}    
    kE2 = {}    
    kE3 = {}    
    kE4 = {}
    
    //constants
    infusion_amount = {}    
    dose_interval = {} // time interval between doses in days    
    dose_end = {} // end of treatment day    
    one_our = 1/24     
    first_dose = {} // time of first dose in days
    const Dpls := infusion_amount;
    
    // events    
    
    end
'''.format(_init_drug_periphery, _init_drug_lung, _init_met_alanine, _init_met_NMP, _init_met_NTP,
           _double_first_dose, _k_p_rate, _k_p_prime_rate, _k01_rate, _k10_rate, _k12_rate, _k23_rate, _k34_rate,
           _kE0_rate,
           _kE1_rate, _kE2_rate, _kE3_rate, _kE4_rate, _infusion_amount, _dose_interval, _eot, _treatment_start)

    drug_dosig_model_vars = ["Dpls", "Dperi", "Dlung", "Mala", "Mnmp", "Mntp"]

    return dosingmodel_str, drug_dosig_model_vars


def set_cell_drug_metabolization(_init_met_alanine, _init_met_NMP, _init_met_NTP,
                                 _k12_rate, _k23_rate, _k34_rate, _kE2_rate, _kE3_rate, _kE4_rate):
    """
    Antimony model string generator for drug metabolization in cells.
    To change parameters do so on the DrugDosingInputs. Parameter descriptions are also in DrugDosingInputs
    :param _init_met_alanine:
    :param _init_met_NMP:
    :param _init_met_NTP:
    :param _k12_rate:
    :param _k23_rate:
    :param _k34_rate:
    :param _kE2_rate:
    :param _kE3_rate:
    :param _kE4_rate:
    :return:
    """

    dosingmodel_str = '''
        model dosingmodel()
        //Time is in days!
    
        //infusion    
        //J0: -> Dpls; switch * infusion_amount / one_our // switch = (0,1) to turn on or off, infusion happens over 1h    
        //flow from plasma     
        //J1: Dpls -> ; kE0 * Dpls // elimination    
        //J2: Dpls -> Dperi ; kp * Dpls // to periphery     
        //J3: Dpls -> Dlung ; k01 * Dpls
            
        // flow from periphery    
        //J4: Dperi -> Dpls ; kpp * Dperi // to plasma     
        
        // Drug reactions / flow in lung    
        //J5: Dlung -> Dpls ; k10 * Dlung    
        //J6: Dlung -> Mala ; k12 * Dlung    
        //J7: Dlung -> ; kE1 * Dlung
        
        // Mala reactions    
        J8: Mala -> Mnmp ; k23 * Mala    
        J9: Mala -> ; kE2 * Mala
        
        //Mnmp reactions    
        J10: Mnmp -> Mntp ; k34 * Mnmp
        J11: Mnmp ->  ; kE3 * Mnmp
        
        // Mntp reaction    
        J12: Mntp -> ; kE4 * Mntp
        
        //parameters
        // initial conditions     
        
        Mala = {}    
        Mnmp = {}    
        Mntp = {}   
        // rates    
        
        k12 = {}    
        k23 = {}    
        k34 = {}        
        kE2 = {}    
        kE3 = {}    
        kE4 = {}

        end
    '''.format(_init_met_alanine, _init_met_NMP, _init_met_NTP, _k12_rate, _k23_rate, _k34_rate,
               _kE2_rate, _kE3_rate, _kE4_rate)

    return dosingmodel_str


def set_simple_pk_full():
    # time units are H!!!
    simple_pk_str = """
                // Created by libAntimony v2.12.0.3
        function Rate_Law_for_Uptake_1(dose, k, duration)
          dose*k/duration;
        end
        
        Rate_Law_for_Uptake_1 is "Rate Law for Uptake_1"
        
        function Function_for_Uptake(Body, Infusion_duration, Remdes_dose_mol, k_in)
          Rate_Law_for_Uptake_1(Remdes_dose_mol, k_in, Infusion_duration)/Body;
        end
        
        Function_for_Uptake is "Function for Uptake"
        
        
        model *New_Model()
        
          // Compartments and Species:
          compartment Body;
          species GS443902 in Body, GS443902_source in Body, GS443902_sink in Body;
        
          // Rate Rules:
          GS443902_AUC' = GS443902;
        
          // Reactions:
          Uptake: GS443902_source => GS443902; Body*Function_for_Uptake(Body, Infusion_duration, Remdes_dose_mol, k_in);
          Clearance: GS443902 => GS443902_sink; Body*k_out*GS443902;
        
          // Events:
          checkTmax: at 0 after GS443902 > GS443902_Cmax: GS443902_Tmax = time;
          checkCmax: at 0 after GS443902 > GS443902_Cmax: GS443902_Cmax = GS443902 + 1e-9;
          checkC24: at 0 after time == 24: GS443902_C24 = GS443902;
        
          // Species initializations:
          GS443902 = 0;
          GS443902_source = Remdes_dose_mol;
          GS443902_sink = 0;
        
          // Compartment initializations:
          Body = 38.4;
        
          // Variable initializations:
          Remdes_dose_mol = Remdes_dose_mg/1000/Remdes_MW;
          Remdes_dose_mol has unit_6;
          GS443902_Cmax = GS443902;
          GS443902_Cmax has unit_4;
          GS443902_Tmax = 0;
          GS443902_Tmax has unit_9;
          GS443902_C24 = 0;
          GS443902_C24 has unit_4;
          k_in = 1;
          k_in has unit_7;
          k_out = ln(2)/Observed_t1_2;
          k_out has unit_7;
          Observed_t1_2 = 30.4;
          Observed_t1_2 has unit_9;
          Remdes_MW = 602.585;
          Remdes_MW has unit_1;
          GS443902_AUC = 0;
          GS443902_AUC has unit_8;
          Remdes_dose_mg = 200;
          Remdes_dose_mg has unit_10;
          Infusion_duration = 1;
          Infusion_duration has unit_9;
        
          // Other declarations:
          var GS443902_Cmax, GS443902_Tmax, GS443902_C24, GS443902_AUC;
          const Body, Remdes_dose_mol, k_in, k_out, Observed_t1_2, Remdes_MW, Remdes_dose_mg;
          const Infusion_duration;
        
          // Unit definitions:
          unit substance = mole;
          unit unit_0 = 1 / 3600e2 second;
          unit unit_1 = gram / mole;
          unit unit_2 = 1 mole * 3600e2 second / litre;
          unit unit_4 = mole / litre;
          unit unit_5 = 1e-3 gram;
          unit unit_6 = mole;
          unit unit_3 = 3600e2 second;
          unit length = metre;
          unit area = metre^2;
          unit volume = litre;
          unit time_unit = 360000e2 second;
          unit unit_7 = 1 / 360000e2 second;
          unit unit_8 = 100 mole * 3600e2 second / litre;
          unit unit_9 = 360000e2 second;
          unit unit_10 = 1e-3 gram;
        
          // Display Names:
          unit_0 is "1/h";
          unit_1 is "g/mol";
          unit_2 is "mol/l*h";
          unit_4 is "mol/l";
          unit_5 is "mg";
          unit_6 is "mol";
          unit_3 is "h";
          time_unit is "time";
          unit_7 is "1/(100*h)";
          unit_8 is "100*mol*h/l";
          unit_9 is "100*h";
          unit_10 is "0.001*g";
        end
        
        New_Model is "New Model_1"


    """

    return simple_pk_str

def set_simple_pk_lung():
    # time units are H!!!
    simple_pk_str = """
        // Created by libAntimony v2.12.0.3
        function Rate_Law_for_Uptake_1(dose, k, duration)
          dose*k/duration;
        end

        Rate_Law_for_Uptake_1 is "Rate Law for Uptake_1"


        model *New_Model()

          // Compartments and Species:
          compartment Body;
          species GS443902 in Body, GS443902_source in Body, GS443902_sink in Body;

          // Rate Rules:
          GS443902_AUC' = GS443902;

          // Reactions:
          Uptake: GS443902_source => GS443902; Rate_Law_for_Uptake_1(Remdes_dose_mol, k_in, Infusion_duration);
          Clearance: GS443902 => GS443902_sink; Body*k_out*GS443902;

          // Events:
          checkTmax: at GS443902 > GS443902_Cmax: GS443902_Tmax = time;
          checkCmax: at GS443902 > GS443902_Cmax: GS443902_Cmax = GS443902 + 1e-8;
          checkC24: at time == 24: GS443902_C24 = GS443902;
          Infusion_0_off: at time > ModelValue_18_0: k_in = 0;
          Infus_1_on: at time > 24: k_in = 0.5;
          Infus_1_off: at time > (24 + ModelValue_18_0): k_in = 0;
          Infus_2_on: at time > 48: k_in = 0.5;
          Infus_2_off: at time > (48 + ModelValue_18_0): k_in = 0;
          Infus_3_on: at time > 72: k_in = 0.5;
          Infus_3_off: at time > (72 + ModelValue_18_0): k_in = 0;
          Infus_4_on: at time > 96: k_in = 0.5;
          Infus_4_off: at time > (96 + ModelValue_18_0): k_in = 0;

          // Species initializations:
          GS443902 = 0;
          GS443902_source = Remdes_dose_mol;
          GS443902_sink = 0;

          // Compartment initializations:
          Body = 38.4;

          // Variable initializations:
          Remdes_dose_mol = Remdes_dose_mg/1000/Remdes_MW;
          Remdes_dose_mol has unit_6;
          GS443902_Cmax = GS443902;
          GS443902_Cmax has unit_4;
          GS443902_Tmax = 0;
          GS443902_Tmax has unit_3;
          GS443902_C24 = 0;
          GS443902_C24 has unit_4;
          ModelValue_18_0 = Infusion_duration;
          k_in = 1;
          k_in has unit_0;
          k_out = ln(2)/Observed_t1_2;
          k_out has unit_0;
          Observed_t1_2 = 30.4;
          Observed_t1_2 has unit_3;
          Remdes_MW = 602.585;
          Remdes_MW has unit_1;
          GS443902_AUC = 0;
          GS443902_AUC has unit_2;
          Remdes_dose_mg = 200;
          Remdes_dose_mg has unit_5;
          Infusion_duration = 1;
          Infusion_duration has unit_3;

          // Other declarations:
          var GS443902_Cmax, GS443902_Tmax, GS443902_C24, k_in, GS443902_AUC;
          const Body, Remdes_dose_mol, ModelValue_18_0, k_out, Observed_t1_2, Remdes_MW;
          const Remdes_dose_mg, Infusion_duration;

          // Unit definitions:
          unit substance = mole;
          unit unit_0 = 1 / 3600e2 second;
          unit unit_1 = gram / mole;
          unit unit_2 = 1 mole * 3600e2 second / litre;
          unit time_unit = 3600e2 second;
          unit unit_4 = mole / litre;
          unit unit_5 = 1e-3 gram;
          unit unit_6 = mole;
          unit length = metre;
          unit area = metre^2;
          unit volume = litre;
          unit unit_3 = 3600e2 second;

          // Display Names:
          unit_0 is "1/h";
          unit_1 is "g/mol";
          unit_2 is "mol/l*h";
          time_unit is "time";
          unit_4 is "mol/l";
          unit_5 is "mg";
          unit_6 is "mol";
          unit_3 is "h";
          ModelValue_18_0 is "Initial for Infusion_duration";
        end

        New_Model is "New Model"

    """

    return simple_pk_str

def set_simple_pk_cell():
    # time units are H!!!
    simple_pk_str = """
        // Created by libAntimony v2.12.0.3
        function Rate_Law_for_Uptake_1(dose, k, duration)
          dose*k/duration;
        end

        Rate_Law_for_Uptake_1 is "Rate Law for Uptake_1"


        model *New_Model()

          // Compartments and Species:
          compartment Body;
          species GS443902 in Body, GS443902_source in Body, GS443902_sink in Body;

          // Rate Rules:
          GS443902_AUC' = GS443902;

          // Reactions:
          //Uptake: GS443902_source => GS443902; Rate_Law_for_Uptake_1(Remdes_dose_mol, k_in, Infusion_duration);
          Clearance: GS443902 => GS443902_sink; Body*k_out*GS443902;

          // Events:
          //checkTmax: at GS443902 > GS443902_Cmax: GS443902_Tmax = time;
          //checkCmax: at GS443902 > GS443902_Cmax: GS443902_Cmax = GS443902 + 1e-8;
          //checkC24: at time == 24: GS443902_C24 = GS443902;
          //Infusion_0_off: at time > ModelValue_18_0: k_in = 0;
          //Infus_1_on: at time > 24: k_in = 0.5;
          //Infus_1_off: at time > (24 + ModelValue_18_0): k_in = 0;
          //Infus_2_on: at time > 48: k_in = 0.5;
          //Infus_2_off: at time > (48 + ModelValue_18_0): k_in = 0;
          //Infus_3_on: at time > 72: k_in = 0.5;
          //Infus_3_off: at time > (72 + ModelValue_18_0): k_in = 0;
          //Infus_4_on: at time > 96: k_in = 0.5;
          //Infus_4_off: at time > (96 + ModelValue_18_0): k_in = 0;

          // Species initializations:
          GS443902 = 0;
          GS443902_source = Remdes_dose_mol;
          GS443902_sink = 0;

          // Compartment initializations:
          Body = 38.4;

          // Variable initializations:
          Remdes_dose_mol = Remdes_dose_mg/1000/Remdes_MW;
          Remdes_dose_mol has unit_6;
          GS443902_Cmax = GS443902;
          GS443902_Cmax has unit_4;
          GS443902_Tmax = 0;
          GS443902_Tmax has unit_3;
          GS443902_C24 = 0;
          GS443902_C24 has unit_4;
          ModelValue_18_0 = Infusion_duration;
          k_in = 1;
          k_in has unit_0;
          k_out = ln(2)/Observed_t1_2;
          k_out has unit_0;
          Observed_t1_2 = 30.4;
          Observed_t1_2 has unit_3;
          Remdes_MW = 602.585;
          Remdes_MW has unit_1;
          GS443902_AUC = 0;
          GS443902_AUC has unit_2;
          Remdes_dose_mg = 200;
          Remdes_dose_mg has unit_5;
          Infusion_duration = 1;
          Infusion_duration has unit_3;

          // Other declarations:
          var GS443902_Cmax, GS443902_Tmax, GS443902_C24, k_in, GS443902_AUC;
          const Body, Remdes_dose_mol, ModelValue_18_0, k_out, Observed_t1_2, Remdes_MW;
          const Remdes_dose_mg, Infusion_duration;

          // Unit definitions:
          unit substance = mole;
          unit unit_0 = 1 / 3600e2 second;
          unit unit_1 = gram / mole;
          unit unit_2 = 1 mole * 3600e2 second / litre;
          unit time_unit = 3600e2 second;
          unit unit_4 = mole / litre;
          unit unit_5 = 1e-3 gram;
          unit unit_6 = mole;
          unit length = metre;
          unit area = metre^2;
          unit volume = litre;
          unit unit_3 = 3600e2 second;

          // Display Names:
          unit_0 is "1/h";
          unit_1 is "g/mol";
          unit_2 is "mol/l*h";
          time_unit is "time";
          unit_4 is "mol/l";
          unit_5 is "mg";
          unit_6 is "mol";
          unit_3 is "h";
          ModelValue_18_0 is "Initial for Infusion_duration";
        end

        New_Model is "New Model"

    """

    return simple_pk_str


class DrugDosingModelSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements drug dosing regimen
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        import Models.DrugDosingModel.DrugDosingInputs as DrugDosingInputs
        BatchRunLib.apply_external_multipliers(__name__, DrugDosingInputs)
        self.drug_dosing_model_key = drug_dosing_model_key

        if not use_simple_pk:
            if constant_drug_concentration:
                self.set_drug_model_string = set_cst_drug_ddm_string
            else:
                self.set_drug_model_string = set_default_ddm_string

            self.set_control_model_string = full_ddm_for_testing
        else:
            self.set_drug_model_string = set_simple_pk_full
            self.set_control_model_string = set_simple_pk_full

        self.plot_ddm_data = plot_ddm_data_freq > 0
        self.write_ddm_data = write_ddm_data_freq > 0

        self.max_avail4 = 2.32417475e-01 * dose  # see comment just before steppable definition

        self.hill_k = active_met_ic50

        self.drug_model_string = None

        self.ddm_vars = None

        self.drug_metabolization_string = None

        self.control_string = None

        self.rmax = None

        self.vr_model_name = ViralInfectionVTMLib.vr_model_name

        self.ddm_rr = None

        self.control_rr = None
        self.active_component = None

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

        if not use_simple_pk:

            self.active_component = 'Mntp'
            # set model string

            self.drug_model_string, self.ddm_vars = self.set_drug_model_string(
                Drug_pls, Drug_peri, Drug_lung, Ala_met, NMP_met, NTP_met,
                first_dose_doubler, kp, kpp, k01, k10, k12, k23, k34, kE0, kE1, kE2, kE3, kE4, dose, dose_interval,
                dose_end,
                first_dose)
            self.control_string, _ = self.set_control_model_string(
                Drug_pls, Drug_peri, Drug_lung, Ala_met, NMP_met, NTP_met,
                first_dose_doubler, kp, kpp, k01, k10, k12, k23, k34, kE0, kE1, kE2, kE3, kE4, dose, dose_interval,
                dose_end,
                first_dose)
            # init sbml
            self.add_free_floating_antimony(model_string=self.drug_model_string, step_size=days_2_mcs,
                                            model_name='drug_dosing_model')
            self.ddm_rr = self.get_roadrunner_for_single_antimony('drug_dosing_model')

            self.add_free_floating_antimony(model_string=self.control_string, step_size=days_2_mcs,
                                            model_name='drug_dosing_control')
            self.control_rr = self.get_roadrunner_for_single_antimony('drug_dosing_control')

            self.drug_metabolization_string = set_cell_drug_metabolization(Ala_met, NMP_met, NTP_met, k12, k23, k34,
                                                                           kE2,
                                                                           kE3, kE4)

            for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                self.add_antimony_to_cell(model_string=self.drug_metabolization_string,
                                          model_name='drug_metabolization',
                                          cell=cell, step_size=days_2_mcs)
        else:
            self.drug_model_string = self.set_simple_pk_full()
            self.active_component = 'GS443902'
            for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                self.add_antimony_to_cell(model_string=self.drug_model_string,
                                          model_name='drug_metabolization',
                                          cell=cell, step_size=hour_2_mcs)





        if prophylactic_treatment:
            # to be able to write the data from prophylaxis I put the prophylactic code in the
            # data steppable. May not be elegant but it works
            # this DOES MEAN that if the write step is not included prophylaxis won't work
            pass

        if sanity_run:
            self.rmax = replicating_rate
        else:
            self.rmax = self.get_rmax(self.sbml.drug_dosing_control['Mntp'])

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
            rate = days_2_mcs * cell.sbml.drug_metabolization['k12']  # k12 is in units of /day!!!!!!!!!!!!!
            u = rate * self.sbml.drug_dosing_model['Dlung'] / len(self.cell_list_by_type(
                self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED))
            total += u
            uptakes.append((cell.id, u))

        return uptakes, total

    def do_prodrug_metabolization(self):
        """

        :return:
        """
        if self.sbml.drug_dosing_model['Dlung'] <= 0:
            # print(self.sbml.drug_dosing_model['Dlung'])
            self.sbml.drug_dosing_model['Dlung'] = 0
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

            if total < self.sbml.drug_dosing_model['Dlung']:
                for cid, u in uptakes:
                    cell = self.fetch_cell_by_id(cid)
                    cell.sbml.drug_metabolization['Mala'] += u
                    # print(cell.id, cell.sbml.drug_metabolization['Available1'])
                    self.sbml.drug_dosing_model['Dlung'] -= u

            else:
                total = self.sbml.drug_dosing_model['Dlung']
                u = self.sbml.drug_dosing_model['Dlung'] / len(self.cell_list_by_type(self.INFECTED,
                                                                                      self.VIRUSRELEASING,
                                                                                      self.UNINFECTED))
                # print('equal uptake = ', u)
                for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                    cell.sbml.drug_metabolization['Mala'] += u
                    self.sbml.drug_dosing_model['Dlung'] -= u

        else:
            # regular cc3d uptake
            secreter = self.get_field_secretor("prodrug")
            total = 0
            for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
                rate = days_2_mcs * cell.sbml.drug_metabolization['k12']
                uptake = secreter.uptakeInsideCellTotalCount(cell, 9e9, rate)  # uptake at rate 'rate' without max
                # value for uptake (9e99)
                cell.sbml.drug_metabolization['Mala'] += abs(uptake.tot_amount)
                self.sbml.drug_dosing_model['Dlung'] -= abs(uptake.tot_amount)
                total += abs(uptake.tot_amount)
        # print('result:', self.sbml.drug_dosing_model['Dlung'], total)

        # print('control:', self.sbml.drug_dosing_control['Dlung'], days_2_mcs * self.sbml.drug_dosing_control['J6'])
        # print('result/control:', self.sbml.drug_dosing_model['Dlung'] / self.sbml.drug_dosing_control['Dlung'],
        #       total / (days_2_mcs * self.sbml.drug_dosing_control['J6']))
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

                cell.dict['rmax'] = self.get_rmax(cell.sbml.drug_metabolization['Mntp'])
                if cell.type != self.UNINFECTED:
                    vr_model = getattr(cell.sbml, self.vr_model_name)
                    vr_model.replicating_rate = cell.dict['rmax']
                    # print(vr_model.replicating_rate)

                time = cell.sbml.drug_metabolization['Time']
            # print('time', time, self.sbml.drug_dosing_model['Time'])

        # for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING):
        #     vr_model = getattr(cell.sbml, self.vr_model_name)
        #     vr_model.replicating_rate = self.rmax

        # ViralInfectionVTMLib.step_sbml_model_cell(cell=cell)
    def simple_pk_step(self, mcs):
        # with the simple pk each cell will have its own pk model running in themselves
        for cell in self.cell_list_by_type(self.INFECTED, self.VIRUSRELEASING, self.UNINFECTED):
            self.timestep_cell_sbml('drug_metabolization', cell)
            cell.dict['rmax'] = self.get_rmax(cell.sbml.drug_metabolization[self.active_component])
            if cell.type != self.UNINFECTED:
                vr_model = getattr(cell.sbml, self.vr_model_name)
                vr_model.replicating_rate = cell.dict['rmax']
        pass

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


class ProdrugDiffusionController(ViralInfectionVTMSteppableBasePy):
    """
    Responsible for controlling prodrug diffusion
    """

    def __init__(self, frequency=1):
        if not diffusing_drug:
            frequency = 9999
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

        self.drug_dosing_model_key = drug_dosing_model_key

        self.step = None
        # self.ddm_rr = None

    def start(self):
        if diffusing_drug:
            self.step = self.diff_step
            self.get_xml_element('prodr_dc').cdata = prodrug_diff_coef_au
            # self.ddm_rr = self.shared_steppable_vars[self.drug_dosing_model_key].ddm_rr
        else:
            self.step = self.empty_step

    def empty_step(self, mcs):
        pass

    def diff_step(self, mcs):

        # 1) compare Dlung in sbml to total amount in simulation (calculate the difference)
        # 1 note) the difference I want is sbml(dlung) - total_rmds_in_sim. that way the sign of the difference is the
        # sign of the flux
        # 2) set cc3d diffusion flux to difference (times rates and whatnot)

        self.get_xml_element('lower_flux').Value = 0.0  # just to be safe

        sim_rmds = self.get_total_prodrug()

        # let me do step 0 first

        dlung = self.sbml.drug_dosing_model['Dlung']

        flux = (dlung - sim_rmds) / (self.dim.x * self.dim.y)  # all pixels get / loose the same share of dlung

        print(sim_rmds, dlung, flux)
        if flux > 0:
            self.get_xml_element('lower_flux').Value = flux
        else:
            if abs(flux) > sim_rmds / (self.dim.x * self.dim.y):
                flux = - sim_rmds / (self.dim.x * self.dim.y)
                self.get_xml_element('lower_flux').Value = flux

    def get_total_prodrug(self):
        try:
            rmds_tot = self.get_field_secretor("prodrug").totalFieldIntegral()
        except AttributeError:  # Pre-v4.2.1 CC3D
            rmds_tot = 0
            for x, y, z in self.every_pixel():
                rmds_tot += self.field.prodrug[x, y, z]
        return rmds_tot


class DrugDosingDataFieldsPlots(ViralInfectionVTMSteppableBasePy):
    """
    Responsible for plots, extra fields and data handling
    """

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)

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

        self.total_virus_released = 0

        self.__flush_counter = 1

        if self.write_ddm_data:
            self.data_files = {'ddm_data': 'ddm_data.dat', 'ddm_rmax_data': 'ddm_rmax_data.dat',
                               'ddm_tot_RNA_data': 'ddm_tot_RNA_data.dat', 'ddm_mean_RNA_data': 'ddm_mean_RNA_data.dat',
                               'ddm_total_viral_production_data': 'ddm_total_viral_production_data.dat'}
            self.ddm_data = {'ddm_data': {}, 'ddm_rmax_data': {}, 'ddm_tot_RNA_data': {}, 'ddm_mean_RNA_data': {},
                             'ddm_total_viral_production_data': {}}

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

        colors = ['blue', 'red', 'green', 'yellow', 'white', 'magenta']
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
                self.shared_steppable_vars['rmax'] = get_rmax(self.mvars, self.sbml.drug_dosing_model['Mala'])
                self.shared_steppable_vars['pre_sim_time'] = number_of_prophylactic_steps - i
                if self.write_ddm_data:
                    self.do_writes(0)
            if self.write_ddm_data:
                self.flush_stored_outputs()
                self.__flush_counter -= 1

    def get_metabolites_in_cell(self, cell):
        # print(cell.sbml.drug_metabolization['Available1'])
        metabolites = [cell.sbml.drug_metabolization[x] for x in self.mvars.ddm_vars[3:]]
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

        self.ddm_data_win.add_data_point(self.mvars.ddm_vars[0], s_to_mcs * mcs / 60 / 60,
                                         self.sbml.drug_dosing_model[self.mvars.ddm_vars[0]])
        self.ddm_data_win.add_data_point(self.mvars.ddm_vars[1], s_to_mcs * mcs / 60 / 60,
                                         self.sbml.drug_dosing_model[self.mvars.ddm_vars[1]])
        self.ddm_data_win.add_data_point(self.mvars.ddm_vars[2], s_to_mcs * mcs / 60 / 60,
                                         self.sbml.drug_dosing_model[self.mvars.ddm_vars[2]])

        total_mets = self.get_total_metabolites_in_cells()

        for i, x in enumerate(self.mvars.ddm_vars[3:]):
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

    def get_ddm_data_list(self):

        d = [self.sbml.drug_dosing_model[x] for x in self.mvars.ddm_vars[:3]]

        total_mets = self.get_total_metabolites_in_cells()
        for m in total_mets:
            d.append(np.sum(m))
        return d

    def do_writes(self, mcs):

        if prophylactic_treatment and mcs == 0:
            time = mcs - self.shared_steppable_vars['pre_sim_time']
        else:
            time = mcs

        if time < 0:
            self.ddm_data['ddm_data'][time] = [self.sbml.drug_dosing_control[x] for x in self.mvars.ddm_vars]

        if time >= 0:
            self.ddm_data['ddm_data'][time] = self.get_ddm_data_list()

            mean, _ = self.get_mean_std_rmax()
            # self.ddm_data['ddm_rmax_data'][time] = [self.shared_steppable_vars['rmax']]
            self.ddm_data['ddm_rmax_data'][time] = [mean]

            rna_list = self.get_rna_array()

            self.ddm_data['ddm_tot_RNA_data'][mcs] = [np.sum(rna_list)]
            self.ddm_data['ddm_mean_RNA_data'][mcs] = [np.mean(rna_list)]

            self.ddm_data['ddm_total_viral_production_data'][mcs] = [self.total_virus_released]

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
        self.total_virus_released += self.shared_steppable_vars['total_virus_release_this_mcs']
        if self.plot_ddm_data and mcs % plot_ddm_data_freq == 0:
            self.do_plots(mcs)
        if self.write_ddm_data and mcs % write_ddm_data_freq == 0:
            self.do_writes(mcs)

    def on_stop(self):
        self.finish()

    def finish(self):
        if self.write_ddm_data:
            self.flush_stored_outputs()
