# coding=utf-8
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

import random

from nCoVToolkit import nCoVUtils
from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase
import ViralInfectionVTMSteppables as MainSteppables
from . import IFNInputs

# Module specific references
IFN_model_name = 'IFN_model'
viral_replication_model_name = 'ifn_viral_replication_model'
ifn_field_name = 'IFNe'
#todo: how to set virus field name
virus_field_name = 'Virus'
# todo: remove simulation control from this module before release
hours_to_simulate = 40.0  # 10 in the original model
ifn_model_vars = ["IFN", "STATP", "IRF7", "IRF7P"]
viral_replication_model_vars = ["H", "V"]

#Steppable keys
ifn_signaling_key = 'ifn_signaling_steppable'
ifn_release_key = 'ifn_release_steppable'  # IFNReleaseSteppable
ifn_field_initializer_key = 'ifn_field_initializer_steppable'  # IFNFieldInitializerSteppable
ifn_sim_data_key = 'ifn_sim_data_steppable'  # SimDataSteppable
ifn_cell_initializer_key = 'ifn_cell_initializer_steppable'
ifn_viral_internalization_key = 'ifn_viral_internalization_steppable'
ifn_ecplise_phase_key = 'ifn_ecplise_phase_steppable'
ifn_viral_release_key = 'ifn_viral_release_steppable'
ifn_viral_death_key = 'ifn_viral_death_steppable'
ifn_virus_field_initializer_key = 'ifn_virus_field_initializer_steppable'
ifn_plaque_assay_key = 'ifn_plaque_assay_steppable'


class IFNSteppableBase(nCoVSteppableBase):
    # todo: basic docstring for IFNSteppableBase
    # todo: document per derived class which of the properties and methods are actually used
    #  e.g., self.register_type corresponds to registration of a cell type for something particular to a derived class

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        # Internal Data
        self._uninfected_type_name = ''
        self._infected_type_name = ''
        self._virus_releasing_type_name = ''
        self._dead_type_name = ''
        self._virus_field_name = ''
        self._ifn_field_name = ''
        self._registered_types = []

    @property
    def uninfected_type_name(self):
        """
        Uninfected cell type name
        """
        return self._uninfected_type_name

    @uninfected_type_name.setter
    def uninfected_type_name(self, _name: str):
        self._uninfected_type_name = _name

    @property
    def infected_type_name(self):
        """
        Infected cell type name
        """
        return self._infected_type_name

    @infected_type_name.setter
    def infected_type_name(self, _name: str):
        self._infected_type_name = _name

    @property
    def virus_releasing_type_name(self):
        """
        Virus-releasing type name
        """
        return self._virus_releasing_type_name

    @virus_releasing_type_name.setter
    def virus_releasing_type_name(self, _name: str):
        self._virus_releasing_type_name = _name

    @property
    def dead_type_name(self):
        """
        Name of the dead cell type
        """
        return self._dead_type_name

    @dead_type_name.setter
    def dead_type_name(self, _name: str):
        self._dead_type_name = _name

    @property
    def virus_field_name(self):
        """
        Name of the virus field
        """
        return self._virus_field_name

    @virus_field_name.setter
    def virus_field_name(self, _name: str):
        self._virus_field_name = _name

    @property
    def ifn_field_name(self):
        """
        Name of the interferon field
        """
        return self._ifn_field_name

    @ifn_field_name.setter
    def ifn_field_name(self, _name: str):
        self._ifn_field_name = _name

    def register_type(self, _name: str):
        if _name not in self._registered_types:
            self._registered_types.append(_name)

    def unregister_type(self, _name: str):
        self._registered_types.remove(_name)

    @property
    def registered_type_ids(self):
        return [getattr(self, x.upper()) for x in self._registered_types]

    @property
    def infected_type_id(self) -> int:
        """
        Id of the infected cell type according to a cc3d simulation
        """
        return getattr(self, self._infected_type_name.upper())

    @property
    def virus_releasing_type_id(self) -> int:
        """
        Id of the virus-releasing cell type according to a cc3d simulation
        """
        return getattr(self, self._virus_releasing_type_name.upper())

    @property
    def dead_type_id(self) -> int:
        """
        Id of the dead cell type according to a cc3d simulation
        """
        return getattr(self, self._dead_type_name.upper())


def IFN_model_string():
    """
    Antimony model string generator for IFN intracellular signaling model adopted from "Multiscale Model
    of RNA Virus Replication and Interferon Responses Reveals Factors Controlling Plaque Growth Dynamics"
    Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego, James A. Glazier and Jason E. Shoemaker
    :return: Antimony model string
    """
    model_string = f'''model {IFN_model_name}()
        //Equations
        E2a: -> IFN         ; H*(k11*RIGI*V+k12*(V^n)/(k13+(V^n))+k14*IRF7P)    ;
        E2b: IFN ->         ; k21*IFN                                           ;
        E4a: -> STATP       ; H*k31*IFNe/(k32+k33*IFNe)                         ;
        E4b: STATP ->       ; t3*STATP                                          ;
        E5a: -> IRF7        ; H*(k41*STATP+k42*IRF7P)                           ;
        E5b: IRF7 ->        ; t4*IRF7                                           ;
        E6a: -> IRF7P       ; H*k51*IRF7                                        ;
        E6b: IRF7P ->       ; t5*IRF7P                                          ;

        //Parameters
        k11 = {IFNInputs.k11}    ;
        k12 = {IFNInputs.k12}    ;
        k13 = {IFNInputs.k13}    ;
        k14 = {IFNInputs.k14}    ;
        k21 = {IFNInputs.k21}    ;
        k31 = {IFNInputs.k31}    ;
        k32 = {IFNInputs.k32}    ;
        k33 = {IFNInputs.k33}    ;
        t3  = {IFNInputs.t3}     ;
        k41 = {IFNInputs.k41}    ;
        k42 = {IFNInputs.k42}    ;
        t4  = {IFNInputs.t4}     ;
        k51 = {IFNInputs.k51}    ;
        t5  = {IFNInputs.t5}     ;
        n   = {IFNInputs.n}      ;
        RIGI = {IFNInputs.RIGI}  ;

        // Inputs
        H    = 0.0     ;
        IFNe = 0.0     ;
        V    = 0.0     ;
    end'''
    return model_string


def viral_replication_model_string():
    """
    Antimony model string generator for viral replication model coupled with the intracellular signaling
    model adopted from "Multiscale Model of RNA Virus Replication and Interferon Responses Reveals Factors
    Controlling Plaque Growth Dynamics" Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego,
    James A. Glazier and Jason E. Shoemaker

    :return: Antimony model string
    """
    model_string = f'''model {viral_replication_model_name}()
        //Equations
        E7a: H ->           ; H*k61*V                     ;
        E8a: -> V           ; H*k71*V/(1.0+k72*IFNe*7E-5) ;
        E8b: V ->           ; k73*V                       ;

        //Parameters
        k61 = {IFNInputs.k61} ;
        k71 = {IFNInputs.k71} ;
        k72 = {IFNInputs.k72} ;
        k73 = {IFNInputs.k73} ;

        //Initial Conditions
        V = 0.0      ;
        H = 1.0      ;

        //Inputs
        IFNe = 0.0   ;
        end
    '''
    return model_string


def get_cell_viral_replication_model(cell):
    """
    Convenience method to get the viral replication model of a cell

    :param cell: a cell
    :return: viral replication model instance
    """
    return getattr(cell.sbml, viral_replication_model_name)


def get_cell_ifn_model(cell):
    """
    Convenience method to get the ifn model of a cell

    :param cell: a cell
    :return: ifn model instance
    """
    return getattr(cell.sbml, IFN_model_name)


class IFNSteppable(IFNSteppableBase):
    """
    Implements IFN signaling model by passing information between viral replication model, IFN signaling model
    and reading the amount of extracellular interferon

    """

    unique_key = ifn_signaling_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        # Internal data
        self._registered_types = []
        self._ifn_params = {k: getattr(IFNInputs, k) for k in ifn_model_vars}

        # Initialize default data
        self.ifn_field_name = ifn_field_name
        self.uninfected_type_name = MainSteppables.uninfected_type_name
        self.infected_type_name = MainSteppables.infected_type_name
        self.virus_releasing_type_name = MainSteppables.virus_releasing_type_name
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name)
        self.register_type(MainSteppables.virus_releasing_type_name)

        self.sbml_options = {'relative': 1e-10, 'absolute': 1e-12}

    def start(self):
        # todo: remove control of simulation steps from this module before release
        hours_to_mcs = self.step_period / 60.0 / 60.0
        self.get_xml_element('simulation_steps').cdata = hours_to_simulate / hours_to_mcs
        self.set_sbml_global_options(self.sbml_options)

        # Load IFN sbml model
        # todo: make model parameters of IFN sbml model settable;
        #  redundant this can be done using the IFNInputs
        def ifn_model_fcn(cell):
            return IFN_model_string()

        self.register_ode_model(model_name=IFN_model_name,
                                model_fcn=ifn_model_fcn,
                                cell_types=self._registered_types,
                                step_size=hours_to_mcs)

        # Load viral replication sbml model
        def vr_model_fcn(cell):
            return viral_replication_model_string()

        self.register_ode_model(model_name=viral_replication_model_name,
                                model_fcn=vr_model_fcn,
                                cell_types=self._registered_types,
                                step_size=hours_to_mcs)

    def step(self, mcs):
        # Connects viral replication model and IFN model and read IFNe field
        secretor = self.get_field_secretor(field_name=self._ifn_field_name)
        for cell in self.cell_list_by_type(*self.registered_type_ids):
            ifn_cell_sbml = get_cell_ifn_model(cell)
            virus_cell_sbml = get_cell_viral_replication_model(cell)
            ifn_cell_sbml['H'] = virus_cell_sbml['H']
            ifn_cell_sbml['V'] = virus_cell_sbml['V']
            ifn_seen = secretor.amountSeenByCell(cell)  # Calculate once, apply twice
            ifn_cell_sbml['IFNe'] = ifn_seen
            virus_cell_sbml['IFNe'] = ifn_seen

        # Step the models for all registered types
        self.timestep_ode_model(model_name=IFN_model_name)
        self.timestep_ode_model(model_name=viral_replication_model_name)


class IFNCellInitializerSteppable(MainSteppables.CellInitializerSteppable, IFNSteppableBase):
    """
    Initalizes cells with the MOI corresponding to the training conditions of the original paper: Weaver JJ,
    Shoemaker JE. Mathematical Modeling of RNA Virus Sensing Pathways Reveals Paracrine Signaling as the Primary Factor
    Regulating Excessive Cytokine Production. Processes. 2020 Jun;8(6):719.
    """

    unique_key = ifn_cell_initializer_key

    def __init__(self, frequency=1):
        super().__init__(frequency=frequency)

        # Initialize default data
        self.virus_releasing_type_name = MainSteppables.virus_releasing_type_name
        self.random_infected_fraction(1.0)

    def start(self):
        super().start()

        for cell in self.cell_list_by_type(self.infected_type_id):
            self.set_cell_type(cell, self.virus_releasing_type_id)


class IFNViralInternalizationSteppable(MainSteppables.ViralInternalizationSteppable, IFNSteppableBase):
    """
    Convenience steppable to rescale infectivity parameter according to cellularization method proposed in:
    Sego TJ, Aponte-Serrano JO, Gianlupi JF, Glazier JA. Generating Agent-Based Multiscale Multicellular
    Spatiotemporal Models from Ordinary Differential Equations of Biological Systems, with Applications
    in Viral Infection. bioRxiv. 2021 Jan 1.
    """

    unique_key = ifn_viral_internalization_key

    def __init__(self, frequency=1):
        MainSteppables.ViralInternalizationSteppable.__init__(self, frequency)
        IFNSteppableBase.__init__(self)

        # Internal data
        self._registered_types = []
        self._initial_amount_virus = 0.0

        # Initialize default data
        self.virus_releasing_type_name = MainSteppables.virus_releasing_type_name
        self.dead_type_name = MainSteppables.dead_type_name
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name)
        self.register_type(MainSteppables.virus_releasing_type_name)
        self.register_type(MainSteppables.dead_type_name)
        self.initial_amount_virus = 6.9e-8

        # To detect if the user specifies an internalization rate, or if we're calculating it with default data
        self._internalization_rate = None

    def step(self, mcs):
        # Calculate internalization rate if necessary
        if self._internalization_rate is None:
            days_to_mcs = self.step_period / 60.0 / 60.0 / 24.0
            initial_number_registered_cells = len(self.cell_list_by_type(*self.registered_type_ids))
            b = IFNInputs.b * initial_number_registered_cells * days_to_mcs
            self.internalization_rate = b

        secretor = self.get_field_secretor(field_name=self._target_field_name)
        for cell in self.cell_list_by_type(self.uninfected_type_id):
            seen_amount = secretor.amountSeenByCell(cell)
            rate = seen_amount * self._internalization_rate
            if random.random() <= nCoVUtils.ul_rate_to_prob(rate):
                self.set_cell_type(cell, self.infected_type_id)

    def on_set_cell_type(self, cell, old_type):
        if cell.type == self.infected_type_id and old_type == self.uninfected_type_id:
            virus_cell_sbml = get_cell_viral_replication_model(cell)
            virus_cell_sbml['V'] = self._initial_amount_virus

    @property
    def initial_amount_virus(self):
        """
        Initial amnount of unitless virus level when internalization occurs
        """
        return self._initial_amount_virus

    @initial_amount_virus.setter
    def initial_amount_virus(self, _val: str):
        self._initial_amount_virus = _val


class IFNEclipsePhaseSteppable(MainSteppables.EclipsePhaseSteppable, IFNSteppableBase):
    """
    Convenience steppable to rescale ecplise phase parameter
    """

    unique_key = ifn_ecplise_phase_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Initialize default data
        days_to_mcs = self.step_period / 60.0 / 60.0 / 24.0
        self.eclipse_phase = 1.0 / (IFNInputs.k * days_to_mcs)


class IFNViralReleaseSteppable(MainSteppables.ViralReleaseSteppable, IFNSteppableBase):
    """
    Performs virion release according to the intracellular virus level

    The rate of virus release per unit of intracellular virus can be set with the attribute ``release_rate``, in
    time units of hours.
    """

    unique_key = ifn_viral_release_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Internal data
        self._virus_level_scaling_factor = 0.0

        # Initialize default data
        self.virus_level_scaling_factor = 1094460.28
        self.release_rate = IFNInputs.k73

    def step(self, mcs):
        secretor = self.get_field_secretor(field_name=self._target_field_name)
        min_dim = min(self.dim.x, self.dim.y, self.dim.z)
        fact = 1.0
        if min_dim < 3:
            fact = float(min_dim)

        hours_to_mcs = self.step_period / 60.0 / 60.0
        for cell in self.cell_list_by_type(self.virus_releasing_type_id):
            virus_cell_sbml = get_cell_viral_replication_model(cell)
            intracellularVirus = virus_cell_sbml['V']
            p = self.release_rate * hours_to_mcs * intracellularVirus * self._virus_level_scaling_factor

            secretor.secreteInsideCell(cell, p * fact / cell.volume)

    @property
    def virus_level_scaling_factor(self):
        """
        Scaling Factor from unitless virus level to PFU/mL
        """
        return self._virus_level_scaling_factor

    @virus_level_scaling_factor.setter
    def virus_level_scaling_factor(self, _val: float):
        self._virus_level_scaling_factor= _val


class IFNViralDeathSteppable(MainSteppables.ViralDeathSteppable, IFNSteppableBase):
    """
    Performs viral death as a function of the cell's viability determined from the cell' internal
    viral replication model.

    The derived class specified the death rate as a per cell rate, rather than a rate by cell type.

    The rate of cell death can be set with the attribute ``viral_death_rate``, in time units of hours.
    """

    unique_key = ifn_viral_death_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Initialize default data
        self.viral_death_rate = IFNInputs.k61

    def step(self, mcs):
        hours_to_mcs = self.step_period / 60.0 / 60.0
        for cell in self.cell_list_by_type(self.virus_releasing_type_id):
            virus_cell_sbml = get_cell_viral_replication_model(cell)
            H = virus_cell_sbml['H']
            V = virus_cell_sbml['V']
            viral_death_rate = self.viral_death_rate * hours_to_mcs * V * (1 - H)
            pr = nCoVUtils.ul_rate_to_prob(viral_death_rate)
            if random.random() <= pr:
                self.set_cell_type(cell, self.dead_type_id)


class IFNReleaseSteppable(IFNSteppableBase):
    """
    Performs IFN release

    If the simulation domain is quasi-2D, then release will be multiplied by the height of the domain.

    The release rate of IFN per unit of intracellular IFN can be set with the attribute ``release_rate``, in
    time units of hours.
    """

    unique_key = ifn_release_key

    def __init__(self, frequency):
        IFNSteppableBase.__init__(self, frequency)

        self.runBeforeMCS = 1

        # Internal Data
        self._registered_types = []
        self._release_rate = 0.0

        # Initialize default data
        self.ifn_field_name = ifn_field_name
        self.uninfected_type_name = MainSteppables.uninfected_type_name
        self.infected_type_name = MainSteppables.infected_type_name
        self.virus_releasing_type_name = MainSteppables.virus_releasing_type_name
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name)
        self.register_type(MainSteppables.virus_releasing_type_name)
        self.release_rate = IFNInputs.k21

    def step(self, mcs):
        secretor = self.get_field_secretor(field_name=self._ifn_field_name)
        min_dim = min(self.dim.x, self.dim.y, self.dim.z)
        fact = 1.0
        if min_dim < 3:
            fact = float(min_dim)

        hours_to_mcs = self.step_period / 60.0 / 60.0
        for cell in self.cell_list_by_type(*self.registered_type_ids):
            ifn_cell_sbml = get_cell_ifn_model(cell)
            intracellularIFN = ifn_cell_sbml['IFN']
            p = self.release_rate * hours_to_mcs * intracellularIFN
            secretor.secreteInsideCell(cell, p * fact / cell.volume)

    @property
    def release_rate(self):
        return self._release_rate

    @release_rate.setter
    def release_rate(self, _val: float):
        if _val < 0:
            raise ValueError("IFN release must be non-negative")
        self._release_rate = _val


class IFNVirusFieldInitializerSteppable(MainSteppables.VirusFieldInitializerSteppable, IFNSteppableBase):
    """
    Initializes virus field data and properties

    By default, requires CC3DML ids "virus_dc" for virus field diffusion coefficient and "virus_decay" for virus field
    decay. These can be set with set_field_data.
    """

    unique_key = ifn_virus_field_initializer_key

    def start(self):
        min_to_mcs = self.step_period / 60.0
        days_to_mcs = min_to_mcs / 60.0 / 24.0
        self.diffusion_coefficient = IFNInputs.virus_diffusion_coefficient * min_to_mcs / self.voxel_length ** 2
        self.decay_coefficient = IFNInputs.c * days_to_mcs
        MainSteppables.VirusFieldInitializerSteppable.start(self)


class IFNFieldInitializerSteppable(IFNSteppableBase):
    """
    Initializes IFN field data and properties

    By default, requires CC3DML ids "ifn_dc" for virus field diffusion coefficient and "ifn_decay" for virus field
    decay
    """

    # todo: make diffusion and decay coefficients settable in IFNFieldInitializerSteppable

    unique_key = ifn_field_initializer_key

    def __init__(self, frequency):
        IFNSteppableBase.__init__(self, frequency)

        # Internal Data
        self._ifn_diffusion_id = 'ifn_dc'
        self._ifn_decay_id = 'ifn_decay'

        # Initialize Defaut Data
        self.ifn_field_name = ifn_field_name

    def start(self):
        min_to_mcs = self.step_period / 60.0
        hours_to_mcs = min_to_mcs / 60.0
        self.diffusion_coefficient = IFNInputs.IFNe_diffusion_coefficient * min_to_mcs / self.voxel_length ** 2
        self.decay_coefficient = IFNInputs.t2 * hours_to_mcs

    @property
    def diffusion_coefficient(self) -> float:
        return self.get_xml_element(self._ifn_diffusion_id).cdata

    @diffusion_coefficient.setter
    def diffusion_coefficient(self, _val: float):
        if _val <= 0:
            raise ValueError("Diffusion coefficient must be positive")
        self.get_xml_element(self._ifn_diffusion_id).cdata = _val

    @property
    def decay_coefficient(self) -> float:
        return self.get_xml_element(self._ifn_decay_id).cdata

    @decay_coefficient.setter
    def decay_coefficient(self, _val: float):
        if _val < 0:
            raise ValueError("Decay coefficient must be non-negative")
        self.get_xml_element(self._ifn_decay_id).cdata = _val

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        if field_name is not None:
            self.ifn_field_name = field_name
        if diffusion is not None:
            self._ifn_diffusion_id = diffusion
        if decay is not None:
            self._ifn_decay_id = decay

    @property
    def field_secretor(self):
        return self.get_field_secretor(self.ifn_field_name)

    @property
    def field_object(self):
        return getattr(self.field, self.ifn_field_name)


class IFNSimDataSteppable(IFNSteppableBase):
    """
    Plots/writes simulation data of interest

    The frequency of plotting cell population data in Player can be set with the attribute plot_pop_data_freq.
    Plotting is disabled when plot_pop_data_freq is set to zero.

    The frequency of writing ell population data  to file can be set with the attribute write_pop_data_freq.
    Data is written to the cc3d output directory with the name 'ifn_data.dat'.
    Writing is disabled when write_ifn_data_freq is set to zero.

    The frequency of plotting interferon model data in Player can be set with the attribute plot_ifn_data_freq.
    Plotting is disabled when plot_ifn_data_freq is set to zero.

    The frequency of writing interferon model data to file can be set with the attribute write_ifn_data_freq.
    Data is written to the cc3d output directory with the name 'ifn_data.dat'.
    Writing is disabled when write_ifn_data_freq is set to zero.

    The frequency of plotting diffusion field data in Player can be set with the attribute plot_med_diff_data_freq.
    Plotting is disabled when plot_med_diff_data_freq is set to zero.

    The frequency of writing diffusion field data to file can be set with the attribute write_med_diff_data_freq.
    Data is written to the cc3d output directory with the name 'med_diff_data.dat'.
    Writing is disabled when write_med_diff_data_freq is set to zero.
    """

    # todo: make all data file paths unique to framework;

    unique_key = ifn_sim_data_key

    def __init__(self, frequency):
        IFNSteppableBase.__init__(self, frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self.pop_data_win = None
        self.pop_data_path = None
        self.pop_data = dict()

        self.plot_pop_data = IFNInputs.plot_pop_data_freq > 0
        self.write_pop_data = IFNInputs.write_pop_data_freq > 0

        self.ifn_data_win = None
        self.ifn_data_path = None
        self.ifn_data = dict()

        self.plot_ifn_data = IFNInputs.plot_ifn_data_freq > 0
        self.write_ifn_data = IFNInputs.write_ifn_data_freq > 0

        self.med_diff_data_win = None
        self.med_diff_data_path = None
        self.med_diff_data = dict()

        self.plot_med_diff_data = IFNInputs.plot_med_diff_data_freq > 0
        self.write_med_diff_data = IFNInputs.write_med_diff_data_freq > 0

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

        # Internal Data
        self._registered_types = []
        self._initial_number_cells = 0

        # Initialize Default Data
        self.uninfected_type_name = MainSteppables.uninfected_type_name
        self.infected_type_name = MainSteppables.infected_type_name
        self.virus_releasing_type_name = MainSteppables.virus_releasing_type_name
        self.dead_type_name = MainSteppables.dead_type_name
        self.virus_field_name = virus_field_name
        self.ifn_field_name = ifn_field_name
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name)
        self.register_type(MainSteppables.virus_releasing_type_name)

    def start(self):
        # Initialize population data plot if requested
        if self.plot_pop_data:
            self.pop_data_win = self.add_new_plot_window(title='Population data',
                                                         x_axis_title='Time (hrs)',
                                                         y_axis_title='Numer of cells',
                                                         x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True,
                                                         config_options={'legend': True})

            self.pop_data_win.add_plot(self._uninfected_type_name, style='Dots', color='blue', size=5)
            self.pop_data_win.add_plot(self._infected_type_name, style='Dots', color='red', size=5)
            self.pop_data_win.add_plot(self._virus_releasing_type_name, style='Dots', color='green', size=5)
            self.pop_data_win.add_plot(self._dead_type_name, style='Dots', color='yellow', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath('pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

                self.pop_data[-1] = ['Time', self._uninfected_type_name, self._infected_type_name,
                                     self._virus_releasing_type_name, self._dead_type_name]

        # Initialize ifn data plot if requested
        if self.plot_ifn_data:
            colors = ['magenta', 'blue', 'green', 'purple']
            for i in range(len(ifn_model_vars)):
                attr_name = 'ifn_data_win' + ifn_model_vars[i]
                new_window = self.add_new_plot_window(title=ifn_model_vars[i],
                                                      x_axis_title='Time (hrs)',
                                                      y_axis_title='Variables', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True,
                                                      config_options={'legend': True})
                new_window.add_plot(ifn_model_vars[i], style='Dots', color=colors[i], size=5)
                setattr(self, attr_name, new_window)

            colors = ['yellow', 'white']
            for i in range(len(viral_replication_model_vars)):
                attr_name = 'ifn_data_win' + viral_replication_model_vars[i]
                new_window = self.add_new_plot_window(title=viral_replication_model_vars[i],
                                                      x_axis_title='Time (hrs)',
                                                      y_axis_title='Variables', x_scale_type='linear',
                                                      y_scale_type='linear',
                                                      grid=True,
                                                      config_options={'legend': True})
                new_window.add_plot(viral_replication_model_vars[i], style='Dots', color=colors[i], size=5)
                setattr(self, attr_name, new_window)

            attr_name = 'ifn_data_win' + self._ifn_field_name
            new_window = self.add_new_plot_window(title=self._ifn_field_name,
                                                  x_axis_title='Time (hrs)',
                                                  y_axis_title='Variables', x_scale_type='linear',
                                                  y_scale_type='linear',
                                                  grid=True,
                                                  config_options={'legend': True})
            new_window.add_plot(self._ifn_field_name, style='Dots', color='red', size=5)
            setattr(self, attr_name, new_window)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_ifn_data:
                self.ifn_data_path = Path(self.output_dir).joinpath('ifn_data.dat')
                with open(self.ifn_data_path, 'w'):
                    pass

                self.ifn_data[-1] = ['Time'] + ifn_model_vars + viral_replication_model_vars + [self._ifn_field_name]

        # Initialize med diff data plot if requested
        if self.plot_med_diff_data:
            self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                              x_axis_title='Time (hrs)',
                                                              y_axis_title='Number of diffusive species per volume',
                                                              x_scale_type='linear',
                                                              y_scale_type='log',
                                                              grid=True,
                                                              config_options={'legend': True})

            self.med_diff_data_win.add_plot(self._virus_field_name, style='Dots', color='red', size=5)
            self.med_diff_data_win.add_plot(self._ifn_field_name, style='Dots', color='green', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_med_diff_data:
                self.med_diff_data_path = Path(self.output_dir).joinpath('med_diff_data.dat')
                with open(self.med_diff_data_path, 'w'):
                    pass

                self.med_diff_data[-1] = ['Time', self._virus_field_name, self._ifn_field_name]

    def step(self, mcs):
        if mcs == 0:
            self.initial_number_cells = len(self.cell_list_by_type(*self.registered_type_ids)) + len(
                self.cell_list_by_type(self.dead_type_id))

        if self.simdata_steppable is None:
            self.simdata_steppable = self.shared_steppable_vars[ifn_sim_data_key]

        plot_pop_data = self.plot_pop_data and mcs % IFNInputs.plot_pop_data_freq == 0
        plot_ifn_data = self.plot_ifn_data and mcs % IFNInputs.plot_ifn_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % IFNInputs.plot_med_diff_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % IFNInputs.write_pop_data_freq == 0
            write_ifn_data = self.write_ifn_data and mcs % IFNInputs.write_ifn_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % IFNInputs.write_med_diff_data_freq == 0
        else:
            write_pop_data = False
            write_ifn_data = False
            write_med_diff_data = False

        hours_to_mcs = self.step_period / 60.0 / 60.0

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.uninfected_type_id))
            num_cells_infected = len(self.cell_list_by_type(self.infected_type_id))
            num_cells_virusreleasing = len(self.cell_list_by_type(self.virus_releasing_type_id))
            num_cells_dying = len(self.cell_list_by_type(self.dead_type_id))

            # Plot population data plot if requested
            if plot_pop_data:
                self.pop_data_win.add_data_point(self._uninfected_type_name,  mcs * hours_to_mcs, num_cells_uninfected)
                self.pop_data_win.add_data_point(self._infected_type_name,  mcs * hours_to_mcs, num_cells_infected)
                self.pop_data_win.add_data_point(self._virus_releasing_type_name,  mcs * hours_to_mcs, num_cells_virusreleasing)
                self.pop_data_win.add_data_point(self._dead_type_name,  mcs * hours_to_mcs, num_cells_dying)

            # Write population data to file if requested
            if write_pop_data:
                self.pop_data[mcs] = [mcs * hours_to_mcs,
                                      num_cells_uninfected,
                                      num_cells_infected,
                                      num_cells_virusreleasing,
                                      num_cells_dying]

        # Plotting and writing average values of IFN model variables
        if plot_ifn_data or write_ifn_data:
            data_dict = [mcs * hours_to_mcs]
            total_vars = {k: 0.0 for k in ifn_model_vars + viral_replication_model_vars}
            L = 0.0

            for cell in self.cell_list_by_type(*self.registered_type_ids):
                L += 1.0

                cell_sbml = get_cell_ifn_model(cell)
                for model_var in ifn_model_vars:
                    total_vars[model_var] += cell_sbml[model_var]

                cell_sbml = get_cell_viral_replication_model(cell)
                for model_var in viral_replication_model_vars:
                    total_vars[model_var] += cell_sbml[model_var]

            for model_var in ifn_model_vars + viral_replication_model_vars:
                mean_var = total_vars[model_var] / L
                data_dict.append(mean_var)
                if plot_ifn_data:
                    window_name = getattr(self, 'ifn_data_win' + model_var)
                    window_name.add_data_point(model_var, mcs * hours_to_mcs, mean_var)

            # Measure amount of IFNe in the Field
            secretor = self.get_field_secretor(field_name=self._ifn_field_name)
            total_var = 0
            for cell in self.cell_list_by_type(*self.registered_type_ids):
                total_var += secretor.amountSeenByCell(cell)
            if plot_ifn_data:
                window_name = getattr(self, 'ifn_data_win' + self._ifn_field_name)
                window_name.add_data_point(self._ifn_field_name, mcs * hours_to_mcs,
                                           total_var / self._initial_number_cells)
            data_dict.append(total_var / self._initial_number_cells)

            if write_ifn_data:
                self.ifn_data[mcs] = data_dict

        if plot_med_diff_data or write_med_diff_data:
            # Gather total diffusive amounts
            try:
                med_viral_total = self.get_field_secretor(field_name=self._virus_field_name).totalFieldIntegral()
                med_ifn_total = self.get_field_secretor(field_name=self._ifn_field_name).totalFieldIntegral()
            except AttributeError:  # Pre-v4.2.1 CC3D
                med_viral_total = 0.0
                med_ifn_total = 0.0
                field_virus = getattr(self.field, self._virus_field_name)
                field_ifn = getattr(self.field, self._ifn_field_name)
                for x, y, z in self.every_pixel():
                    med_viral_total += field_virus[x, y, z]
                    med_ifn_total += field_ifn[x, y, z]

            # Plot total diffusive viral amount if requested
            if plot_med_diff_data:
                if med_viral_total > 0:
                    self.med_diff_data_win.add_data_point(self._virus_field_name, mcs * hours_to_mcs, med_viral_total)
                if med_ifn_total > 0:
                    self.med_diff_data_win.add_data_point(self._ifn_field_name, mcs * hours_to_mcs, med_ifn_total)

            # Write total diffusive viral amount if requested
            if write_med_diff_data:
                self.med_diff_data[mcs] = [mcs * hours_to_mcs, med_viral_total, med_ifn_total]

        self.flush_stored_outputs()

    def on_stop(self):
        self.finish()

    def finish(self):
        self.flush_stored_outputs()

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage
        :return: None
        """
        if self.write_pop_data and self.pop_data:
            with open(self.pop_data_path, 'a') as fout:
                fout.write(MainSteppables.SimDataSteppable.data_output_string(self, self.pop_data))
                self.pop_data.clear()

        if self.write_ifn_data and self.ifn_data:
            with open(self.ifn_data_path, 'a') as fout:
                fout.write(MainSteppables.SimDataSteppable.data_output_string(self, self.ifn_data))
                self.ifn_data.clear()

        if self.write_med_diff_data and self.med_diff_data:
            with open(self.med_diff_data_path, 'a') as fout:
                fout.write(MainSteppables.SimDataSteppable.data_output_string(self, self.med_diff_data))
                self.med_diff_data.clear()

    @property
    def uninfected_type_id(self) -> int:
        """
        Id of the uninfected cell type according to a cc3d simulation
        """
        return getattr(self, self._uninfected_type_name.upper())

    @property
    def initial_number_cells(self):
        """
        Count of initial number of epithelial cells for rescaling purposes
        """
        return self._initial_number_cells

    @initial_number_cells.setter
    def initial_number_cells(self, _val: int):
        self._initial_number_cells = _val

class IFNPlaqueAssaySteppable(IFNSteppableBase):
    """
    Initializes, plots and writes data for plaque assay simulations

    The frequency of plotting plaque assay data in Player can be set with the attribute plot_plaque_assay_data_freq.
    Plotting is disabled when plot__plaque_assay_data_freq is set to zero.

    The frequency of writing ell population data  to file can be set with the attribute write_plaque_assay_data_freq.
    Data is written to the cc3d output directory with the name 'plaque_assay_data.dat'.
    Writing is disabled when write_plaque_assay_data_freq is set to zero.
    """

    unique_key = ifn_plaque_assay_key

    def __init__(self, frequency):
        IFNSteppableBase.__init__(self, frequency)

        # Reference to SimDataSteppable
        self.plaque_assay_data_steppable = None

        self.plaque_assay_data_win = None
        self.plaque_assay_data_path = None
        self.plaque_assay_data = dict()

        self.plot_plaque_assay_data = IFNInputs.plot_plaque_assay_data_freq > 0
        self.write_plaque_assay_data = IFNInputs.write_plaque_assay_data_freq > 0

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

        # Initialize Default Data
        self.infected_type_name = MainSteppables.infected_type_name
        self.virus_releasing_type_name = MainSteppables.virus_releasing_type_name
        self.dead_type_name = MainSteppables.dead_type_name

    def start(self):
        #todo: initialize single infected cell in the center of the simulation domain

        # Initialize population data plot if requested
        if self.plot_plaque_assay_data:
            self.plaque_assay_data_win = self.add_new_plot_window(title='Plaque Data',
                                                         x_axis_title='Time (hrs)',
                                                         y_axis_title='Radial Distance',
                                                         x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=True,
                                                         config_options={'legend': True})

            self.plaque_assay_data_win.add_plot(self._infected_type_name, style='Dots', color='red', size=5)
            self.plaque_assay_data_win.add_plot(self._virus_releasing_type_name, style='Dots', color='green', size=5)
            self.plaque_assay_data_win.add_plot(self._dead_type_name, style='Dots', color='yellow', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_plaque_assay_data:
                self.plaque_assay_data_path = Path(self.output_dir).joinpath('plaque_assay_data.dat')
                with open(self.plaque_assay_data_path, 'w'):
                    pass

                self.plaque_assay_data[-1] = ['Time',self._infected_type_name,self._virus_releasing_type_name,
                                              self._dead_type_name]

    def step(self, mcs):
        if self.plaque_assay_data_steppable is None:
            self.plaque_assay_data_steppable = self.shared_steppable_vars[ifn_plaque_assay_key]

        plot_plaque_assay_data = self.plot_plaque_assay_data and mcs % IFNInputs.plot_plaque_assay_data_freq == 0
        if self.output_dir is not None:
            write_plaque_assay_data = self.write_plaque_assay_data and mcs % IFNInputs.write_plaque_assay_data_freq == 0
        else:
            write_plaque_assay_data = False


        hours_to_mcs = self.step_period / 60.0 / 60.0

        if plot_plaque_assay_data or write_plaque_assay_data:
            import numpy as np
            # Measure area occupied by cells and assume its a circle
            volume_infected = 0.0
            for cell in self.cell_list_by_type(self.infected_type_id):
                volume_infected += cell.volume
            avg_infected_radius = np.sqrt(volume_infected/np.pi)

            volume_virus_releasing = 0.0
            for cell in self.cell_list_by_type(self.virus_releasing_type_id):
                volume_virus_releasing += cell.volume
            avg_virus_releasing_radius = np.sqrt(volume_virus_releasing/np.pi)

            volume_dead = 0.0
            for cell in self.cell_list_by_type(self.dead_type_id):
                volume_dead += cell.volume
            avg_dead_radius = np.sqrt(volume_dead/np.pi)

            # Plot plaque assay data if requested
            if plot_plaque_assay_data:
                self.plaque_assay_data_win.add_data_point(self._infected_type_name, mcs * hours_to_mcs,
                                                 avg_infected_radius)
                self.plaque_assay_data_win.add_data_point(self._virus_releasing_type_name, mcs * hours_to_mcs,
                                                 avg_virus_releasing_radius)
                self.plaque_assay_data_win.add_data_point(self._dead_type_name, mcs * hours_to_mcs, avg_dead_radius)

            # Write population data to file if requested
            if write_plaque_assay_data:
                self.plaque_assay_data[mcs] = [mcs * hours_to_mcs,
                                      avg_infected_radius,
                                      avg_virus_releasing_radius,
                                      avg_dead_radius]
        self.flush_stored_outputs()

    def on_stop(self):
        self.finish()

    def finish(self):
        self.flush_stored_outputs()

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage
        :return: None
        """
        if self.write_plaque_assay_data and self.plaque_assay_data:
            with open(self.plaque_assay_data_path, 'a') as fout:
                fout.write(MainSteppables.SimDataSteppable.data_output_string(self, self.plaque_assay_data))
                self.plaque_assay_data.clear()