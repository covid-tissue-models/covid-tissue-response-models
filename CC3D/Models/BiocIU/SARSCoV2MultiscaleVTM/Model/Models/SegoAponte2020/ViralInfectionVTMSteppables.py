###############################################################################################################
# To cite this model please use the following:
#
# T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
# James M. Osborne, Ellen M. Quardokus, Richard K. Plemper, James A. Glazier,
# "A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and
# immune response in epithelial tissues and its application to drug therapy timing and effectiveness",
# bioRxiv 2020.04.27.064139
###############################################################################################################

import math
import numpy as np

from . import ViralInfectionVTMBasePy
from . import ViralInfectionVTMLib
from . import ViralInfectionVTMModelInputs

# Import main steppables module
import ViralInfectionVTMSteppables as MainSteppables

# Import toolkit
from nCoVToolkit import nCoVUtils

ViralInfectionVTMSteppableBasePy = ViralInfectionVTMBasePy.ViralInfectionVTMSteppableBasePy

# Name of generic immune cell type
immune_type_name = 'Immunecell'
# Name of generic cytokine field
cytokine_field_name = 'cytokine'
# Name of generic oxidative agent field
oxidator_field_name = 'oxidator'

rng = np.random  # alias for random number generators (rng)


class CellsInitializerSteppable(ViralInfectionVTMSteppableBasePy, MainSteppables.CellInitializerSteppable):
    """
    Initializes epithelial sheet with a single infected cell

    Requires ViralReplicationSteppable registered beforehand to automatically instantiate cells with viral
    replication model. Otherwise, viral replication model needs to be instantiated per cell after initialization
    """

    unique_key = ViralInfectionVTMLib.cell_initializer_steppable_key

    def __init__(self, frequency=1):
        ViralInfectionVTMSteppableBasePy.__init__(self, frequency)
        MainSteppables.CellInitializerSteppable.__init__(self, frequency)

        # Initialize default data
        self.single_infected_cell()

    def start(self):
        ViralInfectionVTMSteppableBasePy.start(self)
        MainSteppables.CellInitializerSteppable.start(self)

        for cell in self.cell_list_by_type(self.uninfected_type_id, self.infected_type_id):
            cell.targetVolume = ViralInfectionVTMModelInputs.cell_volume
            cell.lambdaVolume = ViralInfectionVTMModelInputs.volume_lm
            if ViralInfectionVTMLib.vrl_key in cell.dict.keys():
                ViralInfectionVTMLib.pack_viral_replication_variables(cell)

        # Set initially infected cell data
        for cell in self.cell_list_by_type(self.infected_type_id):
            var_unpacking = ViralInfectionVTMLib.vr_cell_dict_to_sym[ViralInfectionVTMLib.vrm_unpacking]
            getattr(cell.sbml, self.vr_model_name)[var_unpacking] = 1.0


class VirusFieldInitializerSteppable(MainSteppables.VirusFieldInitializerSteppable):

    unique_key = ViralInfectionVTMLib.virus_field_initializer_key

    def __init__(self, frequency=1):
        super().__init__(frequency=frequency)

    def start(self):

        self.diffusion_coefficient = ViralInfectionVTMModelInputs.virus_dc
        self.decay_coefficient = ViralInfectionVTMModelInputs.virus_decay


class ViralReplicationSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral replication module

    Assigns a callback to instantiation of new epithelial cells to assign a new intracellular model

    Additional cell types can be registered and unregisterd with register_type and unregister_type, respectively.
    Cell types that are registered with the keyword argument `infected` set to True will have their intracellular model
    integrated in time every simulation step. All registrations must occur before start is called.

    By default, module uninfected, infected and virus-releasing types are registered, with the latter two registered as
    infected.

    The criterion for the transition from module types infected to virus releasing and virus releasing to dead are
    described by cell_releases and cell_dies, respectively, which can be set/overridden as desired.
    """

    unique_key = ViralInfectionVTMLib.vrm_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self._registered_types = []
        self._infected_types = []

        # Initialize default data
        self.sbml_options = {'relative': 1e-10, 'absolute': 1e-12}
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name, infected=True)
        self.register_type(MainSteppables.virus_releasing_type_name, infected=True)

    def start(self):
        # Load model
        self.set_sbml_global_options(self.sbml_options)

        self.register_ode_model(model_name=self.vr_model_name,
                                model_fcn=self.epithelial_model_fcn_generator(),
                                cell_types=self._registered_types,
                                step_size=ViralInfectionVTMModelInputs.vr_step_size)

    def step(self, mcs):
        if self.simdata_steppable is None:
            try:
                self.simdata_steppable: SimDataSteppable = self.shared_steppable_vars[SimDataSteppable.unique_key]
            except KeyError:
                pass

        # Do viral model
        for cell in self.cell_list_by_type(*self.infected_type_ids):
            # Step the model for this cell
            ViralInfectionVTMLib.step_sbml_model_cell(cell=cell)
            # Pack state variables into cell dictionary
            ViralInfectionVTMLib.pack_viral_replication_variables(cell=cell)

            # Test for infected -> virus-releasing
            if self.cell_releases(cell):
                self.set_cell_type(cell, self.virus_releasing_type_id)

            # Test for virus-releasing -> dead
            if self.cell_dies(cell):
                self.kill_cell(cell=cell)
                if self.simdata_steppable is not None:
                    self.simdata_steppable.track_death_viral()

    def cell_releases(self, _cell):
        return _cell.type == self.infected_type_id and \
               _cell.dict[ViralInfectionVTMLib.vrm_assembled] > ViralInfectionVTMModelInputs.cell_infection_threshold

    def cell_dies(self, _cell):
        return _cell.type == self.virus_releasing_type_id and \
               np.random.random() < nCoVUtils.hill_equation(_cell.dict[ViralInfectionVTMLib.vrm_assembled],
                                                            ViralInfectionVTMModelInputs.diss_coeff_uptake_apo,
                                                            ViralInfectionVTMModelInputs.hill_coeff_uptake_apo)

    def epithelial_model_fcn_generator(self):
        def model_fcn(_cell):
            _cell.dict[ViralInfectionVTMLib.vrl_key] = True
            for k in [ViralInfectionVTMLib.vrm_unpacking,
                      ViralInfectionVTMLib.vrm_replicating,
                      ViralInfectionVTMLib.vrm_packing,
                      ViralInfectionVTMLib.vrm_assembled,
                      ViralInfectionVTMLib.vrm_uptake]:
                if k not in _cell.dict.keys():
                    _cell.dict[k] = 0.0
            return self.viral_replication_model_string(ViralInfectionVTMModelInputs.unpacking_rate,
                                                       ViralInfectionVTMModelInputs.replicating_rate,
                                                       ViralInfectionVTMModelInputs.r_half,
                                                       ViralInfectionVTMModelInputs.translating_rate,
                                                       ViralInfectionVTMModelInputs.packing_rate,
                                                       ViralInfectionVTMModelInputs.secretion_rate,
                                                       _cell.dict[ViralInfectionVTMLib.vrm_unpacking],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_replicating],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_packing],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_assembled],
                                                       _cell.dict[ViralInfectionVTMLib.vrm_uptake])
        return model_fcn

    def register_type(self, _type_name: str, infected: bool = False):
        self._registered_types.append(_type_name)
        if infected:
            self._infected_types.append(_type_name)

    def unregister_type(self, _type_name: str):
        self._registered_types.remove(_type_name)
        if _type_name in self._infected_types:
            self._infected_types.remove(_type_name)

    @property
    def infected_type_ids(self):
        return [getattr(self, x.upper()) for x in self._infected_types]

    @property
    def registered_type_ids(self):
        return [getattr(self, x.upper()) for x in self._registered_types]

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback. Ensures that necessary data is available and properly set
        """
        if _new_cell.type in self.registered_type_ids and ViralInfectionVTMLib.vrl_key not in _new_cell.dict.keys():
            _new_cell.dict[ViralInfectionVTMLib.vrl_key] = False


class ViralInternalizationSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral internalization module
    """
    # todo: make initial_unbound_receptors settable on ViralInternalizationSteppable

    unique_key = ViralInfectionVTMLib.vim_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Initialize default data
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)

    def do_cell_internalization(self, cell, viral_amount_com):
        """
        Describes internalization criterion

        :param cell: a cell
        :param viral_amount_com: total amount of virus associated with a cell
        :return: boolean of whether internalization occurs, and corresponding uptake amount
        :rtype: (bool, float)
        """
        if cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] == 0:
            return False, 0.0

        _k = ViralInfectionVTMModelInputs.kon * cell.volume / ViralInfectionVTMModelInputs.koff
        init_unbound_recept = ViralInfectionVTMModelInputs.initial_unbound_receptors
        hill_coeff_uptake_pr = ViralInfectionVTMModelInputs.hill_coeff_uptake_pr
        receptors = cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key]
        diss_coeff_uptake_pr = (init_unbound_recept / 2.0 / _k / receptors) ** (1.0 / hill_coeff_uptake_pr)
        uptake_probability = nCoVUtils.hill_equation(viral_amount_com,
                                                     diss_coeff_uptake_pr,
                                                     hill_coeff_uptake_pr)

        cell_does_uptake = np.random.rand() < uptake_probability
        rate_coeff_uptake_pr = ViralInfectionVTMModelInputs.rate_coeff_uptake_pr
        uptake_amount = ViralInfectionVTMModelInputs.s_to_mcs / rate_coeff_uptake_pr * uptake_probability

        if cell_does_uptake and cell.type == self.uninfected_type_id:
            self.set_cell_type(cell, self.infected_type_id)

        return cell_does_uptake, uptake_amount

    def update_cell_receptors(self, cell, receptors_increment):
        """
        Updates the number of unbound receptors on a cell's surface

        :param cell: a cell
        :param receptors_increment: increment to number of receptors
        :return: None
        """
        cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] = max(
            cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] + receptors_increment, 0.0)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback. Ensures that necessary data is available and properly set
        """
        if _new_cell.type == self.uninfected_type_id and \
                ViralInfectionVTMLib.unbound_receptors_cellg_key not in _new_cell.dict.keys():
            _new_cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key] = \
                ViralInfectionVTMModelInputs.initial_unbound_receptors


class ViralSecretionSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements viral release module

    This steppable uses a callback per registered cell type that is called on each cell of the registered type every
    time step to calculate interactions with the module virus field.

    Every secretion callback has the signature secr_func(_steppable, _cell, _mcs)
    - _steppable (ViralSecretionSteppable): self
    - _cell (cc3d.cpp.CompuCell.CellG): a cell
    - _mcs (int): the current simulation step

    Callbacks can be registered for a cell type with register_secretor_by_type, and unregistered for a cell type with
    unregister_secretor_by_type.

    By default, module uninfected, infected and virus-releasing types are assigned the callback secr_func_epithelial.

    By default, requires a diffusion field with module virus field name. If all types registered with the callback
    secr_func_epithelial are unregistered, then this requirement is no longer relevant.

    By default, requires ViralInternalizationSteppable for determining/implementing internalization events. This can be
    customized by setting `vim_steppable` with any object that implements do_cell_internalization and
    update_cell_receptors with the same signature as those defined by ViralInternalizationSteppable. If all types
    registered with the callback secr_func_epithelial are unregistered, then this requirement is no longer relevant.
    """

    unique_key = ViralInfectionVTMLib.vrs_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        # Reference to ViralInternalizationSteppable
        self.vim_steppable = None

        self._secretors_by_type = {}

        # Initialize default data
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_virus_field_name(MainSteppables.virus_field_name)
        self.register_secretor_by_type(MainSteppables.uninfected_type_name, self.secr_func_epithelial)
        self.register_secretor_by_type(MainSteppables.infected_type_name, self.secr_func_epithelial)
        self.register_secretor_by_type(MainSteppables.virus_releasing_type_name, self.secr_func_epithelial)

    def start(self):
        if ViralInfectionVTMModelInputs.track_model_variables:
            for attribute_name in [ViralInfectionVTMLib.vrm_unpacking,
                                   ViralInfectionVTMLib.vrm_replicating,
                                   ViralInfectionVTMLib.vrm_packing,
                                   ViralInfectionVTMLib.vrm_assembled,
                                   ViralInfectionVTMLib.vrm_secretion]:
                field_name = attribute_name.replace(ViralInfectionVTMLib.module_prefix, '')
                self.track_cell_level_scalar_attribute(field_name=field_name, attribute_name=attribute_name)

    def step(self, mcs):
        if self.vim_steppable is None:
            try:
                self.vim_steppable: ViralInternalizationSteppable = \
                    self.shared_steppable_vars[ViralInternalizationSteppable.unique_key]
            except KeyError:
                pass

        for type_name, secr_func in self._secretors_by_type.items():
            for cell in self.cell_list_by_type(getattr(self, type_name.upper())):
                secr_func(self, cell, mcs)

    @staticmethod
    def secr_func_epithelial(self, _cell, _mcs):
        # Evaluate probability of cell uptake of viral particles from environment
        # If cell isn't infected, it changes type to infected here if uptake occurs
        secretor = self.get_field_secretor(self._virus_field_name)
        viral_amount_com = self.total_seen_field(field=self.virus_field, cell=_cell, estimate=True)
        cell_does_uptake, uptake_amount = self.vim_steppable.do_cell_internalization(_cell, viral_amount_com)
        if cell_does_uptake:
            uptake = secretor.uptakeInsideCellTotalCount(_cell, 1E12, uptake_amount / _cell.volume)
            _cell.dict[ViralInfectionVTMLib.vrm_uptake] = abs(uptake.tot_amount)
            self.vim_steppable.update_cell_receptors(
                cell=_cell,
                receptors_increment=-_cell.dict[ViralInfectionVTMLib.vrm_uptake]*ViralInfectionVTMModelInputs.s_to_mcs)
            ViralInfectionVTMLib.set_viral_replication_cell_uptake(cell=_cell,
                                                                   uptake=_cell.dict[ViralInfectionVTMLib.vrm_uptake])

        if _cell.type == self.virus_releasing_type_id:
            sec_amount = ViralInfectionVTMLib.get_viral_replication_cell_secretion(cell=_cell)
            secretor.secreteInsideCellTotalCount(_cell, sec_amount / _cell.volume)

    def register_secretor_by_type(self, _type_name: str, _secr_func):
        """
        Register a secretion callback by cell type

        :param _type_name: name of cell type
        :param _secr_func: secretion callback
        :return: None
        """
        self._secretors_by_type[_type_name] = _secr_func

    def unregister_secretor_by_type(self, _type_name: str):
        """
        Unregister a secretion callback by cell type

        :param _type_name: name of cell type
        :return: callback of registered cell type
        """
        return self._secretors_by_type.pop(_type_name)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback. If a cell's type changes to the registered virus-releasing type, then
        release of virus in enabled
        """
        # todo: make secretion_rate settable on ViralSecretionSteppable
        if cell.type == self.virus_releasing_type_id:
            ViralInfectionVTMLib.enable_viral_secretion(cell=cell,
                                                        secretion_rate=ViralInfectionVTMModelInputs.secretion_rate)


class ImmuneCellKillingSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell direct cytotoxicity and bystander effect module

    Cell types that can be killed by contact-mediated interactions with the module immune cell type can be
    registered and unregistered with append_contact_target and remove_contact_target, respectively.

    Cell types that can be killed by the bystander effect can be registered and unregistered with add_bystander_target
    and remove_bystander_target, respectively.

    By default, module infected and virus-releasing types are registered for contact-mediated killing.

    By default, module uninfected, infected and virus-releasing types are registered for bystander effect killing.
    """

    unique_key = ViralInfectionVTMLib.immune_killing_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self._contact_targets = []
        self._bystander_targets = {}

        # Initialize defaults
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.set_immune_type_name(immune_type_name)
        # Infected and virus-releasing types are contact targets
        self.append_contact_target(MainSteppables.infected_type_name)
        self.append_contact_target(MainSteppables.virus_releasing_type_name)
        # Uninfected, infected and virus-releasing types are bystander targets
        self.add_bystander_target(MainSteppables.uninfected_type_name, ViralInfectionVTMModelInputs.bystander_effect)
        self.add_bystander_target(MainSteppables.infected_type_name, ViralInfectionVTMModelInputs.bystander_effect)
        self.add_bystander_target(MainSteppables.virus_releasing_type_name,
                                  ViralInfectionVTMModelInputs.bystander_effect)

    def step(self, mcs):
        if self.simdata_steppable is None:
            try:
                self.simdata_steppable: SimDataSteppable = self.shared_steppable_vars[SimDataSteppable.unique_key]
            except KeyError:
                pass

        killed_cells = []
        for cell in self.cell_list_by_type(*self.contact_target_ids):
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type == self.immune_type_id:
                        self.kill_cell(cell=cell)
                        killed_cells.append(cell)
                        if self.simdata_steppable is not None:
                            self.simdata_steppable.track_death_contact()

        # Bystander Effect
        for cell in killed_cells:
            for neighbor, common_surface_area in self.get_cell_neighbor_data_list(cell):
                if neighbor:
                    if neighbor.type in self.bystander_target_ids:
                        if np.random.random() < self._bystander_targets[self.get_type_name_by_cell(neighbor)]:
                            self.kill_cell(cell=neighbor)
                            if self.simdata_steppable is not None:
                                self.simdata_steppable.track_death_bystander()

    def append_contact_target(self, _name: str):
        """
        Add a cell type target for contact-mediated killing

        :param _name: name of cell type
        :return: None
        """
        self._contact_targets.append(_name)

    def remove_contact_target(self, _name: str):
        """
        Remove a cell type target for contact-mediated killing

        :param _name: name of cell type
        :return: None
        """
        self._contact_targets.remove(_name)

    def contact_target_names(self) -> list:
        return self._contact_targets.copy()

    @property
    def contact_target_ids(self):
        return [getattr(self, x.upper()) for x in self.contact_target_names()]

    def add_bystander_target(self, _name: str, _prob: float):
        """
        Add a cell type target for bystander effect killing

        :param _name: name of cell type
        :param _prob: probability of killing by the bystander effect
        :return: None
        """
        if not 0 <= _prob <= 1:
            raise ValueError('Bystander effect probability must be in [0, 1]')
        self._bystander_targets[_name] = _prob

    def remove_bystander_target(self, _name: str):
        """
        Remove a cell type target for bystander effect killing

        :param _name: name of cell type
        :return: None
        """
        self._bystander_targets.pop(_name)

    @property
    def bystander_targets(self):
        return [x for x in self._bystander_targets.keys()]

    @property
    def bystander_target_ids(self):
        return [getattr(self, x.upper()) for x in self.bystander_targets]


class ChemotaxisSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell chemotaxis module

    By default, immune cells chemotaxis according to the module cytokine field. The name of the field can be set with
    set_target_field_name.
    """

    unique_key = ViralInfectionVTMLib.chemotaxis_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self._target_field_name = ''
        self._immune_type_name = ''
        self._lamda_chemotaxis = 0.0

        # Initialize defaults
        self.set_immune_type_name(immune_type_name)
        self.set_target_field_name(cytokine_field_name)
        self.set_lamda_chemotaxis(ViralInfectionVTMModelInputs.lamda_chemotaxis)

    def start(self):
        for cell in self.cell_list_by_type(self.immune_type_id):
            self.add_cell_chemotaxis_data(cell=cell)

    def step(self, mcs):
        field = self.target_field
        for cell in self.cell_list_by_type(self.immune_type_id):

            cd = self.chemotaxisPlugin.getChemotaxisData(cell, self._target_field_name)
            if cell.dict[ViralInfectionVTMLib.activated_cellg_key]:
                cd.setLambda(self._lamda_chemotaxis / (1.0 + field[cell.xCOM, cell.yCOM, 1]))
            else:
                cd.setLambda(0)

    def set_target_field_name(self, _name: str):
        self._target_field_name = _name

    def set_lamda_chemotaxis(self, _val: float):
        self._lamda_chemotaxis = _val

    @property
    def target_field(self):
        return getattr(self.field, self._target_field_name)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback.

        If a cell is of registered immune cell type, then its chemotaxis data is initialized
        """
        if _new_cell.type == self.immune_type_id:
            self.add_cell_chemotaxis_data(cell=_new_cell)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback.

        If a cell is of registered immune cell type, then its chemotaxis data is initialized
        """
        if cell.type == self.immune_type_id:
            self.add_cell_chemotaxis_data(cell=cell)

    def add_cell_chemotaxis_data(self, cell):
        cd = self.chemotaxisPlugin.addChemotaxisData(cell, self._target_field_name)
        cd.assignChemotactTowardsVectorTypes([self.MEDIUM])


class ImmuneCellSeedingSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell seeding and removal of immune cell recruitment module

    Seeded immune cells are initialized with the model parameters of this module in their dictionary

    The name of the field by which immune cells are seeded can be set with set_seeding_field_name

    The name of the field by which immune cells chemotaxis can be set with set_chemotaxis_field_name

    By default, ImmuneRecruitmentSteppable to determine seeding and removal probabilities.
    If customizing, this functionality can be replaced by setting the attribute ir_steppable with any object that
    provides methods get_immune_seeding_prob() and get_immune_removal_prob() that return the probability of seeding and
    removing an immune cell, respectively
    If deriving, this functionality can be replaced by overriding get_immune_seeding_prob and get_immune_removal_prob
    to describe the tests for seeding and removing an immune cell, respectively
    """

    unique_key = ViralInfectionVTMLib.immune_seeding_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to ImmuneResponseSteppable
        self.ir_steppable = None

        self._seeding_field_name = ''
        self._chemotaxis_field_name = ''

        # Initialize defaults
        self.set_immune_type_name(immune_type_name)
        self.set_seeding_field_name(MainSteppables.virus_field_name)
        self.set_chemotaxis_field_name(cytokine_field_name)

    def start(self):
        cell_diameter = ViralInfectionVTMModelInputs.cell_diameter

        for iteration in range(int(ViralInfectionVTMModelInputs.initial_immune_seeding)):
            cell = True
            while cell:
                xi = np.random.randint(0, self.dim.x - 2 * cell_diameter)
                yi = np.random.randint(0, self.dim.y - 2 * cell_diameter)
                for x in range(xi, xi + int(cell_diameter)):
                    for y in range(yi, yi + int(cell_diameter)):
                        cell = self.cell_field[x, y, 1]
                        break
                cell = False
            cell = self.new_immune_cell_in_time(ck_production=ViralInfectionVTMModelInputs.max_ck_secrete_im,
                                                ck_consumption=ViralInfectionVTMModelInputs.max_ck_consume)
            self.cell_field[x:x + int(cell_diameter), y:y + int(cell_diameter), 1] = cell
            cell.targetVolume = ViralInfectionVTMModelInputs.cell_volume
            cell.lambdaVolume = ViralInfectionVTMModelInputs.volume_lm

    def step(self, mcs):
        if self.ir_steppable is None:
            try:
                self.ir_steppable: ImmuneRecruitmentSteppable = \
                    self.shared_steppable_vars[ImmuneRecruitmentSteppable.unique_key]
            except KeyError:
                pass

        for cell in self.cell_list_by_type(self.immune_type_id):
            p_immune_dying = np.random.random()
            if p_immune_dying < self.get_immune_removal_prob():
                cell.targetVolume = 0.0

        p_immune_seeding = np.random.random()
        if p_immune_seeding < self.get_immune_seeding_prob():
            open_space = True
            viral_concentration = 0
            cell_diameter = ViralInfectionVTMModelInputs.cell_diameter
            for iteration in range(10):
                radius = 10
                length = 0
                while length <= radius:
                    xi = np.random.randint(0, self.dim.x - 2 * ViralInfectionVTMModelInputs.cell_diameter)
                    yi = np.random.randint(0, self.dim.y - 2 * ViralInfectionVTMModelInputs.cell_diameter)
                    length = np.sqrt((self.dim.x // 2 - xi) ** 2 + (self.dim.y // 2 - yi) ** 2)
                for x in range(xi, xi + int(ViralInfectionVTMModelInputs.cell_diameter)):
                    for y in range(yi, yi + int(ViralInfectionVTMModelInputs.cell_diameter)):
                        cell = self.cell_field[x, y, 1]
                        if cell:
                            open_space = False
                            break
                if open_space:
                    concentration_iteration = self.seeding_field[xi, yi, 1]
                    if concentration_iteration >= viral_concentration:
                        viral_concentration = concentration_iteration
                        x_seed = xi
                        y_seed = yi
            if open_space:
                cell = self.new_immune_cell_in_time(ck_production=ViralInfectionVTMModelInputs.max_ck_secrete_im,
                                                    ck_consumption=ViralInfectionVTMModelInputs.max_ck_consume)

                self.cell_field[x_seed:x_seed + int(cell_diameter), y_seed:y_seed + int(cell_diameter), 1] = cell

                cell.targetVolume = ViralInfectionVTMModelInputs.cell_volume
                cell.lambdaVolume = ViralInfectionVTMModelInputs.volume_lm

    def set_seeding_field_name(self, _name: str):
        """
        Set the name of the seeding field

        :param _name: seeding field name
        :return: None
        """
        self._seeding_field_name = _name

    @property
    def seeding_field(self):
        return getattr(self.field, self._seeding_field_name)

    def set_chemotaxis_field_name(self, _name: str):
        """
        Set the name of the target chemotaxis field

        :param _name: chemotaxis field name
        :return: None
        """
        self._chemotaxis_field_name = _name

    def get_immune_seeding_prob(self):
        """
        Get the current probability of seeding an immune cell

        :return: seeding probability
        """
        if self.ir_steppable is not None:
            return self.ir_steppable.get_immune_seeding_prob()
        return 0.0

    def get_immune_removal_prob(self):
        """
        Get the current probability of removing an immune cell

        :return: removal probability
        """
        if self.ir_steppable is not None:
            return self.ir_steppable.get_immune_removal_prob()
        return 0.0


class SimDataSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Plots/writes simulation data of interest
    """

    unique_key = ViralInfectionVTMLib.simdata_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.vrm_data_win = None
        self.vrm_data_path = None
        self.vrm_data = dict()
        # The viral replication model of this cell is tracked and plotted/recorded
        self.vrm_tracked_cell = None

        self.vim_data_win = None
        self.vim_data_path = None
        self.vim_data = dict()

        self.pop_data_win = None
        self.pop_data_path = None
        self.pop_data = dict()

        self.med_diff_data_win = None
        self.med_diff_data_path = None
        self.med_diff_data = dict()

        self.ir_data_win = None
        self.ir_data_path = None
        self.ir_data = dict()

        self.spat_data_win = None
        self.spat_data_path = None
        self.spat_data = dict()

        self.death_data_win = None
        self.death_data_path = None
        self.death_data = dict()

        self.plot_vrm_data = ViralInfectionVTMModelInputs.plot_vrm_data_freq > 0
        self.write_vrm_data = ViralInfectionVTMModelInputs.write_vrm_data_freq > 0

        self.plot_vim_data = ViralInfectionVTMModelInputs.plot_vim_data_freq > 0
        self.write_vim_data = ViralInfectionVTMModelInputs.write_vim_data_freq > 0

        self.plot_pop_data = ViralInfectionVTMModelInputs.plot_pop_data_freq > 0
        self.write_pop_data = ViralInfectionVTMModelInputs.write_pop_data_freq > 0

        self.plot_med_diff_data = ViralInfectionVTMModelInputs.plot_med_diff_data_freq > 0
        self.write_med_diff_data = ViralInfectionVTMModelInputs.write_med_diff_data_freq > 0
        self.med_diff_key = "MedDiff"

        self.plot_ir_data = ViralInfectionVTMModelInputs.plot_ir_data_freq > 0
        self.write_ir_data = ViralInfectionVTMModelInputs.write_ir_data_freq > 0
        self.ir_key = "ImmuneResp"
        self.ir_steppable = None

        self.plot_spat_data = ViralInfectionVTMModelInputs.plot_spat_data_freq > 0
        self.write_spat_data = ViralInfectionVTMModelInputs.write_spat_data_freq > 0

        self.plot_death_data = ViralInfectionVTMModelInputs.plot_death_data_freq > 0
        self.write_death_data = ViralInfectionVTMModelInputs.write_death_data_freq > 0

        # Origin of infection point; if more than one cell is first detected, then measure the mean COM
        # If first infection is far from center of domain, then measurements of infection front probably won't
        # be very useful
        self.init_infect_pt = None

        # Cell death mechanism tracking
        self._death_mech = {'viral': 0,
                            'oxi': 0,
                            'contact': 0,
                            'bystander': 0}

        # For flushing outputs every quarter simulation length
        self.__flush_counter = 1

        # Initialize defaults
        self.set_uninfected_type_name(MainSteppables.uninfected_type_name)
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.set_immune_type_name(immune_type_name)
        self.set_virus_field_name(MainSteppables.virus_field_name)
        self.set_cytokine_field_name(cytokine_field_name)
        self.set_oxidator_field_name(oxidator_field_name)

    def start(self):
        if self.plot_vrm_data:
            self.vrm_data_win = self.add_new_plot_window(title='VRM',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Variables', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=False,
                                                         config_options={'legend': True})

            self.vrm_data_win.add_plot("U", style='Dots', color='blue', size=5)
            self.vrm_data_win.add_plot("R", style='Dots', color='orange', size=5)
            self.vrm_data_win.add_plot("P", style='Dots', color='green', size=5)
            self.vrm_data_win.add_plot("A", style='Dots', color='red', size=5)
            self.vrm_data_win.add_plot("Uptake", style='Dots', color='yellow', size=5)
            self.vrm_data_win.add_plot("Secretion", style='Dots', color='white', size=5)

        if self.plot_vim_data:
            self.vim_data_win = self.add_new_plot_window(title='VIM',
                                                         x_axis_title='MonteCarlo Step (MCS)',
                                                         y_axis_title='Variables', x_scale_type='linear',
                                                         y_scale_type='linear',
                                                         grid=False,
                                                         config_options={'legend': True})

            self.vim_data_win.add_plot("R", style='Dots', color='orange', size=5)

        # Initialize population data plot if requested
        if self.plot_pop_data:
            self.pop_data_win = self.add_new_plot_window(title='Population data',
                                                         x_axis_title='MCS',
                                                         y_axis_title='Numer of cells',
                                                         x_scale_type='linear',
                                                         y_scale_type='log',
                                                         grid=True,
                                                         config_options={'legend': True})

            self.pop_data_win.add_plot("Uninfected", style='Dots', color='blue', size=5)
            self.pop_data_win.add_plot("Infected", style='Dots', color='red', size=5)
            self.pop_data_win.add_plot("VirusReleasing", style='Dots', color='green', size=5)
            self.pop_data_win.add_plot("Dying", style='Dots', color='yellow', size=5)
            self.pop_data_win.add_plot("ImmuneCell", style='Dots', color='white', size=5)
            self.pop_data_win.add_plot("ImmuneCellActivated", style='Dots', color='purple', size=5)

        if self.plot_med_diff_data:
            self.med_diff_data_win = self.add_new_plot_window(title='Total diffusive species',
                                                              x_axis_title='MCS',
                                                              y_axis_title='Number of diffusive species per volume',
                                                              x_scale_type='linear',
                                                              y_scale_type='log',
                                                              grid=True,
                                                              config_options={'legend': True})

            self.med_diff_data_win.add_plot("MedViral", style='Dots', color='red', size=5)
            self.med_diff_data_win.add_plot("MedCyt", style='Dots', color='blue', size=5)
            self.med_diff_data_win.add_plot("MedOxi", style='Dots', color='green', size=5)

        if self.plot_ir_data:
            self.ir_data_win = self.add_new_plot_window(title='Immune Response Model',
                                                        x_axis_title='MCS',
                                                        y_axis_title='State variable S',
                                                        x_scale_type='linear',
                                                        y_scale_type='linear',
                                                        grid=True)

            self.ir_data_win.add_plot(self.ir_key, style='Dots', color='red', size=5)

        if self.plot_spat_data:
            self.spat_data_win = self.add_new_plot_window(title='Spatial data',
                                                          x_axis_title='MCS',
                                                          y_axis_title='',
                                                          x_scale_type='linear',
                                                          y_scale_type='linear',
                                                          grid=True,
                                                          config_options={'legend': True})

            self.spat_data_win.add_plot("DeathComp", style='Dots', color='red', size=5)
            self.spat_data_win.add_plot("InfectDist", style='Dots', color='blue', size=5)

        if self.plot_death_data:
            self.death_data_win = self.add_new_plot_window(title='Death data',
                                                           x_axis_title='MCS',
                                                           y_axis_title='Numer of cells',
                                                           x_scale_type='linear',
                                                           y_scale_type='log',
                                                           grid=True,
                                                           config_options={'legend': True})

            self.death_data_win.add_plot("Viral", style='Dots', color='blue', size=5)
            self.death_data_win.add_plot("OxiField", style='Dots', color='red', size=5)
            self.death_data_win.add_plot("Contact", style='Dots', color='green', size=5)
            self.death_data_win.add_plot("Bystander", style='Dots', color='yellow', size=5)

        # Check that output directory is available
        if self.output_dir is not None:
            from pathlib import Path
            if self.write_vrm_data:
                self.vrm_data_path = Path(self.output_dir).joinpath(ViralInfectionVTMLib.module_prefix + 'vrm_data.dat')
                with open(self.vrm_data_path, 'w'):
                    pass

            if self.write_vim_data:
                self.vim_data_path = Path(self.output_dir).joinpath(ViralInfectionVTMLib.module_prefix + 'vim_data.dat')
                with open(self.vim_data_path, 'w'):
                    pass

            if self.write_pop_data:
                self.pop_data_path = Path(self.output_dir).joinpath(ViralInfectionVTMLib.module_prefix + 'pop_data.dat')
                with open(self.pop_data_path, 'w'):
                    pass

            if self.write_med_diff_data:
                self.med_diff_data_path = Path(self.output_dir).joinpath(
                    ViralInfectionVTMLib.module_prefix + 'med_diff_data.dat')
                with open(self.med_diff_data_path, 'w'):
                    pass

            if self.write_ir_data:
                self.ir_data_path = Path(self.output_dir).joinpath(ViralInfectionVTMLib.module_prefix + 'ir_data.dat')
                with open(self.ir_data_path, 'w'):
                    pass

            if self.write_spat_data:
                self.spat_data_path = Path(self.output_dir).joinpath(
                    ViralInfectionVTMLib.module_prefix + 'spat_data.dat')
                with open(self.spat_data_path, 'w'):
                    pass

            if self.write_death_data:
                self.death_data_path = Path(self.output_dir).joinpath(
                    ViralInfectionVTMLib.module_prefix + 'death_data.dat')
                with open(self.death_data_path, 'w'):
                    pass

    def step(self, mcs):

        plot_pop_data = self.plot_pop_data and mcs % ViralInfectionVTMModelInputs.plot_pop_data_freq == 0
        plot_med_diff_data = self.plot_med_diff_data and mcs % ViralInfectionVTMModelInputs.plot_med_diff_data_freq == 0
        plot_ir_data = self.plot_ir_data and mcs % ViralInfectionVTMModelInputs.plot_ir_data_freq == 0
        plot_vrm_data = self.plot_vrm_data and mcs % ViralInfectionVTMModelInputs.plot_vrm_data_freq == 0
        plot_vim_data = self.plot_vim_data and mcs % ViralInfectionVTMModelInputs.plot_vim_data_freq == 0
        plot_spat_data = self.plot_spat_data and mcs % ViralInfectionVTMModelInputs.plot_spat_data_freq == 0
        plot_death_data = self.plot_death_data and mcs % ViralInfectionVTMModelInputs.plot_death_data_freq == 0
        if self.output_dir is not None:
            write_pop_data = self.write_pop_data and mcs % ViralInfectionVTMModelInputs.write_pop_data_freq == 0
            write_med_diff_data = self.write_med_diff_data and mcs % ViralInfectionVTMModelInputs.write_med_diff_data_freq == 0
            write_ir_data = self.write_ir_data and mcs % ViralInfectionVTMModelInputs.write_ir_data_freq == 0
            write_vrm_data = self.write_vrm_data and mcs % ViralInfectionVTMModelInputs.write_vrm_data_freq == 0
            write_vim_data = self.write_vim_data and mcs % ViralInfectionVTMModelInputs.write_vim_data_freq == 0
            write_spat_data = self.write_spat_data and mcs % ViralInfectionVTMModelInputs.write_spat_data_freq == 0
            write_death_data = self.write_death_data and mcs % ViralInfectionVTMModelInputs.write_death_data_freq == 0
        else:
            write_pop_data = False
            write_med_diff_data = False
            write_ir_data = False
            write_vrm_data = False
            write_vim_data = False
            write_spat_data = False
            write_death_data = False

        if self.vrm_tracked_cell is not None and (plot_vrm_data or write_vrm_data):
            if plot_vrm_data:
                self.vrm_data_win.add_data_point("U", mcs,
                                                 self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_unpacking])
                self.vrm_data_win.add_data_point("R", mcs,
                                                 self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_replicating])
                self.vrm_data_win.add_data_point("P", mcs,
                                                 self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_packing])
                self.vrm_data_win.add_data_point("A", mcs,
                                                 self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_assembled])
                self.vrm_data_win.add_data_point("Secretion", mcs,
                                                 self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_secretion])

            if write_vrm_data:
                self.vrm_data[mcs] = [self.vrm_tracked_cell.id,
                                      self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_unpacking],
                                      self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_replicating],
                                      self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_packing],
                                      self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_assembled],
                                      self.vrm_tracked_cell.dict[ViralInfectionVTMLib.vrm_secretion]]

        if self.vrm_tracked_cell is not None and (plot_vim_data or write_vim_data):
            if plot_vim_data:
                self.vim_data_win.add_data_point(
                    "R", mcs, self.vrm_tracked_cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key])

            if write_vim_data:
                self.vim_data[mcs] = [self.vrm_tracked_cell.id,
                                      self.vrm_tracked_cell.dict[ViralInfectionVTMLib.unbound_receptors_cellg_key]]

        if plot_pop_data or write_pop_data:

            # Gather population data
            num_cells_uninfected = len(self.cell_list_by_type(self.uninfected_type_id))
            num_cells_infected = len(self.cell_list_by_type(self.infected_type_id))
            num_cells_virusreleasing = len(self.cell_list_by_type(self.virus_releasing_type_id))
            num_cells_dying = len(self.cell_list_by_type(self.dead_type_id))
            num_cells_immune = len(self.cell_list_by_type(self.immune_type_id))
            num_cells_immune_act = len([c for c in self.cell_list_by_type(self.immune_type_id)
                                        if c.dict[ViralInfectionVTMLib.activated_cellg_key]])

            # Plot population data plot if requested
            if plot_pop_data:
                if num_cells_uninfected > 0:
                    self.pop_data_win.add_data_point('Uninfected', mcs, num_cells_uninfected)
                if num_cells_infected > 0:
                    self.pop_data_win.add_data_point('Infected', mcs, num_cells_infected)
                if num_cells_virusreleasing > 0:
                    self.pop_data_win.add_data_point('VirusReleasing', mcs, num_cells_virusreleasing)
                if num_cells_dying > 0:
                    self.pop_data_win.add_data_point('Dying', mcs, num_cells_dying)
                if num_cells_immune > 0:
                    self.pop_data_win.add_data_point('ImmuneCell', mcs, num_cells_immune)
                if num_cells_immune_act > 0:
                    self.pop_data_win.add_data_point('ImmuneCellActivated', mcs, num_cells_immune_act)

            # Write population data to file if requested
            if write_pop_data:
                self.pop_data[mcs] = [num_cells_uninfected,
                                      num_cells_infected,
                                      num_cells_virusreleasing,
                                      num_cells_dying,
                                      num_cells_immune,
                                      num_cells_immune_act]

        if plot_med_diff_data or write_med_diff_data:

            # Gather total diffusive amounts
            try:
                med_viral_total = self.get_field_secretor(self._virus_field_name).totalFieldIntegral()
                med_cyt_total = self.get_field_secretor(self._cytokine_field_name).totalFieldIntegral()
                med_oxi_total = self.get_field_secretor(self._oxidator_field_name).totalFieldIntegral()
            except AttributeError:  # Pre-v4.2.1 CC3D
                med_viral_total = 0.0
                med_cyt_total = 0.0
                med_oxi_total = 0.0
                for x, y, z in self.every_pixel():
                    med_viral_total += self.virus_field[x, y, z]
                    med_cyt_total += self.cytokine_field[x, y, z]
                    med_oxi_total += self.oxidator_field[x, y, z]

            # Plot total diffusive viral amount if requested
            if plot_med_diff_data:
                if med_viral_total > 0:
                    self.med_diff_data_win.add_data_point("MedViral", mcs, med_viral_total)
                if med_cyt_total > 0:
                    self.med_diff_data_win.add_data_point("MedCyt", mcs, med_cyt_total)
                if med_oxi_total > 0:
                    self.med_diff_data_win.add_data_point("MedOxi", mcs, med_oxi_total)

            # Write total diffusive viral amount if requested
            if write_med_diff_data:
                self.med_diff_data[mcs] = [med_viral_total,
                                           med_cyt_total,
                                           med_oxi_total]

        if plot_ir_data or write_ir_data:
            if self.ir_steppable is None:
                self.ir_steppable: ImmuneRecruitmentSteppable = self.shared_steppable_vars[
                    ImmuneRecruitmentSteppable.unique_key]

            s_val = self.ir_steppable.get_state_variable_val()

            # Plot state variable S if requested
            if plot_ir_data:
                self.ir_data_win.add_data_point(self.ir_key, mcs, s_val)

            # Write state variable S if requested
            if write_ir_data:
                self.ir_data[mcs] = [s_val]

        if plot_spat_data or write_spat_data:
            # Calculate compactness of dead cell area as total surface area of intefaces between dying and non-dying
            # types in epithelial sheet divided by total volume of dying types
            dead_srf = 0
            dead_vol = 0
            dying_cell_list = self.cell_list_by_type(self.dead_type_id)
            if not dying_cell_list:
                dead_comp = 0
            else:
                for cell in dying_cell_list:
                    dead_vol += cell.volume
                    for neighbor, common_srf in self.get_cell_neighbor_data_list(cell):
                        if neighbor is not None and neighbor.type in [self.uninfected_type_id,
                                                                      self.infected_type_id,
                                                                      self.virus_releasing_type_id]:
                            dead_srf += common_srf

                dead_comp = dead_srf / dead_vol

            # Calculate infection front: max. distance from initial point of infection to all infected cells
            # If no infected cells, distance is -1
            max_infect_dist = -1
            if self.init_infect_pt is None:
                infected_cell_list = self.cell_list_by_type(self.infected_type_id, self.virus_releasing_type_id)
                num_cells_infected = len(infected_cell_list)
                if num_cells_infected > 0:
                    self.init_infect_pt = [0, 0, 0]
                    for cell in infected_cell_list:
                        self.init_infect_pt[0] += cell.xCOM
                        self.init_infect_pt[1] += cell.yCOM

                    self.init_infect_pt[0] /= num_cells_infected
                    self.init_infect_pt[1] /= num_cells_infected

            if self.init_infect_pt is not None:
                for cell in self.cell_list_by_type(self.infected_type_id, self.virus_releasing_type_id):
                    dx = cell.xCOM - self.init_infect_pt[0]
                    dy = cell.yCOM - self.init_infect_pt[1]
                    max_infect_dist = max(max_infect_dist, math.sqrt(dx * dx + dy * dy))

            # Plot spatial data if requested
            #   Infection distance is normalized by average lattice dimension
            if plot_spat_data:
                self.spat_data_win.add_data_point("DeathComp", mcs, dead_comp)
                if max_infect_dist > 0:
                    max_infect_dist_norm = max_infect_dist / ((self.dim.x + self.dim.y) / 2.0)
                    self.spat_data_win.add_data_point("InfectDist", mcs, max_infect_dist_norm)

            # Write spatial data if requested
            if write_spat_data:
                self.spat_data[mcs] = [dead_comp,
                                       max_infect_dist]

        if plot_death_data or write_death_data:
            num_viral = self._death_mech['viral']
            num_oxi = self._death_mech['oxi']
            num_contact = self._death_mech['contact']
            num_bystander = self._death_mech['bystander']

            # Plot death data if requested
            if plot_death_data:
                if num_viral > 0:
                    self.death_data_win.add_data_point("Viral", mcs, num_viral)
                if num_oxi > 0:
                    self.death_data_win.add_data_point("OxiField", mcs, num_oxi)
                if num_contact > 0:
                    self.death_data_win.add_data_point("Contact", mcs, num_contact)
                if num_bystander > 0:
                    self.death_data_win.add_data_point("Bystander", mcs, num_bystander)

            # Write death data if requested
            if write_death_data:
                self.death_data[mcs] = [num_viral,
                                        num_oxi,
                                        num_contact,
                                        num_bystander]

        # Flush outputs at quarter simulation lengths
        if mcs >= int(self.simulator.getNumSteps() / 4 * self.__flush_counter):
            self.flush_stored_outputs()
            self.__flush_counter += 1

    def on_stop(self):
        self.finish()

    def finish(self):
        self.flush_stored_outputs()

    def data_output_string(self, _data: dict):
        """
        Generate string for data output to file from data dictionary

        :param _data: data dictionary; keys are steps, values are lists of data
        :return: output string to write to file
        """
        mcs_list = list(_data.keys())
        mcs_list.sort()
        f_str = ''
        for mcs in mcs_list:
            f_str += f'{mcs}'
            for v in _data[mcs]:
                f_str += f', {v}'
            f_str += '\n'
        return f_str

    def flush_stored_outputs(self):
        """
        Write stored outputs to file and clear output storage

        :return: None
        """
        # Each tuple contains the necessary information for writing a set of data to file
        #   1. Boolean for whether we're writing to file at all
        #   2. The path to write the data to
        #   3. The data to write
        output_info = [(self.write_vrm_data, self.vrm_data_path, self.vrm_data),
                       (self.write_vim_data, self.vim_data_path, self.vim_data),
                       (self.write_pop_data, self.pop_data_path, self.pop_data),
                       (self.write_med_diff_data, self.med_diff_data_path, self.med_diff_data),
                       (self.write_ir_data, self.ir_data_path, self.ir_data),
                       (self.write_death_data, self.death_data_path, self.death_data)]
        for write_data, data_path, data in output_info:
            if write_data:
                with open(data_path, 'a') as fout:
                    fout.write(self.data_output_string(data))
                    data.clear()

    def set_vrm_tracked_cell(self, cell):
        self.vrm_tracked_cell = cell

    def track_death_viral(self):
        self._death_mech['viral'] += 1

    def track_death_oxi_field(self):
        self._death_mech['oxi'] += 1

    def track_death_contact(self):
        self._death_mech['contact'] += 1

    def track_death_bystander(self):
        self._death_mech['bystander'] += 1


class CytokineProductionAbsorptionSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements cytokine production/secretion and immune cell activation module, and sets cytokine field data

    This steppable uses a callback per registered cell type that is called on each cell of the registered type every
    time step to calculate interactions with the module cytokine field.

    Every secretion callback has the signature val = secr_func(_steppable, _cell, _mcs)
    - _steppable (CytokineProductionAbsorptionSteppable): self
    - _cell (cc3d.cpp.CompuCell.CellG): a cell
    - _mcs (int): the current simulation step
    - val (float): total produced cytokine by the cell

    Callbacks can be registered for a cell type with register_secretor_by_type, and unregistered for a cell type with
    unregister_secretor_by_type.

    By default, requires a diffusion field with module cytokine field name and CC3DML ids 'cytokine_dc' and
    'cytokine_decay' for the field diffusion and decay constants, respectively. These can be set with set_field_data.

    By default, module infected and virus-releasing types are assigned the callback secr_func_infected and module
    immune cell type is assigned the callback secr_func_immune.
    """

    unique_key = ViralInfectionVTMLib.cytokine_secretion_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        if ViralInfectionVTMModelInputs.track_model_variables:
            attribute_name = ViralInfectionVTMLib.activated_cellg_key
            field_name = attribute_name.replace(ViralInfectionVTMLib.module_prefix, '')
            self.track_cell_level_scalar_attribute(field_name=field_name, attribute_name=attribute_name)

        self._cytokine_diffusion_id = ''
        self._cytokine_decay_id = ''
        self._secretors_by_type = {}

        # Reference to ImmuneResponseSteppable
        self.ir_steppable = None

        # Initialize defaults
        self.set_infected_type_name(MainSteppables.infected_type_name)
        self.set_virus_releasing_type_name(MainSteppables.virus_releasing_type_name)
        self.set_immune_type_name(immune_type_name)
        self.set_field_data(field_name=cytokine_field_name, diffusion='cytokine_dc', decay='cytokine_decay')
        self.register_secretor_by_type(MainSteppables.infected_type_name, self.secr_func_infected)
        self.register_secretor_by_type(MainSteppables.virus_releasing_type_name, self.secr_func_infected)
        self.register_secretor_by_type(immune_type_name, self.secr_func_immune)

    def start(self):
        # cytokine diff parameters
        self.get_xml_element(self._cytokine_diffusion_id).cdata = ViralInfectionVTMModelInputs.cytokine_dc
        self.get_xml_element(self._cytokine_decay_id).cdata = ViralInfectionVTMModelInputs.cytokine_field_decay

        for cell in self.cell_list_by_type(self.immune_type_id):
            # cytokine production/uptake parameters for immune cells

            cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = ViralInfectionVTMModelInputs.max_ck_secrete_im
            cell.dict[ViralInfectionVTMLib.ck_consumption_cellg_key] = ViralInfectionVTMModelInputs.max_ck_consume

        for cell in self.cell_list_by_type(self.infected_type_id, self.virus_releasing_type_id):
            cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = ViralInfectionVTMModelInputs.max_ck_secrete_infect

    def step(self, mcs):
        if self.ir_steppable is None:
            self.ir_steppable: ImmuneRecruitmentSteppable = \
                self.shared_steppable_vars[ImmuneRecruitmentSteppable.unique_key]

        # Track the total amount added and subtracted to the cytokine field
        total_ck_inc = 0.0

        for type_name, secr_func in self._secretors_by_type.items():
            for cell in self.cell_list_by_type(getattr(self, type_name.upper())):
                total_ck_inc += secr_func(self, cell, mcs)

        self.ir_steppable.increment_total_cytokine_count(total_ck_inc)

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        if field_name is not None:
            self.set_cytokine_field_name(field_name)
        if diffusion is not None:
            self._cytokine_diffusion_id = diffusion
        if decay is not None:
            self._cytokine_decay_id = decay

    @staticmethod
    def secr_func_infected(self, _cell, _mcs):
        ck_secretor = self.get_field_secretor(self._cytokine_field_name)
        viral_load = ViralInfectionVTMLib.get_assembled_viral_load_inside_cell(
            _cell, ViralInfectionVTMModelInputs.vr_step_size)
        produced = _cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] * nCoVUtils.hill_equation(
            viral_load, ViralInfectionVTMModelInputs.ec50_infecte_ck_prod, 2)
        res = ck_secretor.secreteInsideCellTotalCount(_cell, produced / _cell.volume)
        return res.tot_amount

    @staticmethod
    def secr_func_immune(self, _cell, _mcs):
        total_ck_inc = 0.0
        ck_secretor = self.get_field_secretor(self._cytokine_field_name)

        up_res = ck_secretor.uptakeInsideCellTotalCount(
            _cell, _cell.dict[ViralInfectionVTMLib.ck_consumption_cellg_key] / _cell.volume, 0.1)
        # decay seen ck
        _cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] *= ViralInfectionVTMModelInputs.ck_memory_immune

        # uptake ck

        # from POV of secretion uptake is negative
        _cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] -= up_res.tot_amount
        total_ck_inc += up_res.tot_amount
        p_activate = nCoVUtils.hill_equation(_cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key],
                                             ViralInfectionVTMModelInputs.EC50_ck_immune,
                                             2)

        minimum_activated_time = ViralInfectionVTMModelInputs.minimum_activated_time
        if rng.uniform() < p_activate and not _cell.dict[ViralInfectionVTMLib.activated_cellg_key]:

            _cell.dict[ViralInfectionVTMLib.activated_cellg_key] = True
            _cell.dict[ViralInfectionVTMLib.time_activation_cellg_key] = _mcs
        elif (_cell.dict[ViralInfectionVTMLib.activated_cellg_key]
              and _mcs - _cell.dict[ViralInfectionVTMLib.time_activation_cellg_key] > minimum_activated_time):
            _cell.dict[ViralInfectionVTMLib.activated_cellg_key] = False
            _cell.dict[ViralInfectionVTMLib.time_activation_cellg_key] = - 99

        if _cell.dict[ViralInfectionVTMLib.activated_cellg_key]:
            seen_field = self.total_seen_field(self.field.cytokine, _cell)
            produced = _cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] * nCoVUtils.hill_equation(
                seen_field, ViralInfectionVTMModelInputs.ec50_immune_ck_prod, 1)
            sec_res = ck_secretor.secreteInsideCellTotalCount(_cell, produced / _cell.volume)

            total_ck_inc += sec_res.tot_amount

        return total_ck_inc

    def register_secretor_by_type(self, _type_name: str, _secr_func):
        """
        Register a secretion callback by cell type

        :param _type_name: name of cell type
        :param _secr_func: secretion callback
        :return: None
        """
        self._secretors_by_type[_type_name] = _secr_func

    def unregister_secretor_by_type(self, _type_name: str):
        """
        Unregister a secretion callback by cell type

        :param _type_name: name of cell type
        :return: callback of registered cell type
        """
        return self._secretors_by_type.pop(_type_name)

    def on_new_cell(self, _new_cell):
        """
        Implementation of callback.

        If a new cell is of registered infected, virus-releasing or immune types, then its cytokine production
        parameter is set.

        If a new cell is of registered immune cell type, total cytokine uptake variable is initialized as zero and
        activated flag is initialized as False as necessary
        """
        if _new_cell.type == self.immune_type_id:
            if ViralInfectionVTMLib.tot_ck_upt_cellg_key not in _new_cell.dict.keys():
                _new_cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] = 0
            if ViralInfectionVTMLib.activated_cellg_key not in _new_cell.dict.keys():
                _new_cell.dict[ViralInfectionVTMLib.activated_cellg_key] = False
        self.on_set_cell_type(cell=_new_cell, old_type=None)

    def on_set_cell_type(self, cell, old_type):
        """
        Implementation of callback. If a cell's type changes to the registered virus-releasing type, then
        the cytokine production model parameter is updated
        """
        if cell.type in [self.infected_type_id, self.virus_releasing_type_id]:
            # update cytokine param
            cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = ViralInfectionVTMModelInputs.max_ck_secrete_infect
        elif cell.type == self.immune_type_id:
            cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = ViralInfectionVTMModelInputs.max_ck_secrete_im
            cell.dict[ViralInfectionVTMLib.ck_consumption_cellg_key] = ViralInfectionVTMModelInputs.max_ck_consume


class ImmuneRecruitmentSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell recruitment module

    Note that total cytokine is currently tracked elsewhere by counting uptake and secretion, and by applying the
    field decay rate applied to the cytokine field. This is only relevant for periodic and zero-flux boundary conditions
    for the diffusive cytokine field.

    By default, looks for a CC3DML tag 'cytokine_decay' to get the field decay coefficient of cytokine. The decay
    coefficient can instead be set with set_cytokine_decay

    To manage the total cytokine manually, use increment_total_cytokine_count
    """
    # todo: make ImmuneRecruitmentSteppable ODE model step size settable

    unique_key = ViralInfectionVTMLib.ir_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Reference to solver
        self.__rr = None

        # Running value of total cytokine; to be updated externally through accessor
        self.__total_cytokine = 0.0

        self.__ck_decay = None

        # Initialize defaults
        self.set_immune_type_name(immune_type_name)

    def start(self):
        if self.__ck_decay is None:
            self.__ck_decay = float(self.get_xml_element('cytokine_decay').cdata)

        # Initialize model
        self.__init_fresh_recruitment_model()

    def step(self, mcs):

        # Update total count of immune cells
        num_immune_cells = len(self.cell_list_by_type(self.immune_type_id))

        # Apply consumption / transmission decay to running total
        total_cytokine_decayed = self.__total_cytokine * self.__ck_decay
        self.__total_cytokine -= total_cytokine_decayed

        # Update model
        total_cytokine_transmitted = ViralInfectionVTMModelInputs.ir_transmission_coeff * total_cytokine_decayed
        self.update_running_recruitment_model(num_immune_cells, total_cytokine_transmitted)

    def finish(self):
        pass

    def __init_fresh_recruitment_model(self):
        # Generate solver instance
        def model_fcn():
            return ViralInfectionVTMLib.immune_recruitment_model_string(
                ViralInfectionVTMModelInputs.ir_add_coeff,
                ViralInfectionVTMModelInputs.ir_subtract_coeff,
                ViralInfectionVTMModelInputs.ir_delay_coeff,
                ViralInfectionVTMModelInputs.ir_decay_coeff)

        self.register_ode_model(model_name=ViralInfectionVTMLib.ir_model_name,
                                model_fcn=model_fcn,
                                step_size=ViralInfectionVTMModelInputs.vr_step_size)

        # Get reference to solver
        from cc3d.CompuCellSetup import persistent_globals as pg
        for model_name, rr in pg.free_floating_sbml_simulators.items():
            if model_name == ViralInfectionVTMLib.ir_model_name:
                self.__rr = rr

    def update_running_recruitment_model(self, num_immune_cells, total_cytokine):
        self.__rr['numImmuneCells'] = num_immune_cells
        self.__rr['totalCytokine'] = total_cytokine
        self.timestep_ode_model(model_name=ViralInfectionVTMLib.ir_model_name)

    def get_state_variable_val(self):
        return self.__rr['S']

    def get_immune_seeding_prob(self):
        """
        Returns probability of immune cell seeding due to local and global recruitment

        Probability is only non-zero if the state variable *S* is positive, in which case
        the probability is an error function of *S*

        :return: probability of immune cell seeding due to local and global recruitment
        """
        s_val = self.get_state_variable_val()
        if s_val < 0:
            return 0.0
        else:
            return math.erf(ViralInfectionVTMModelInputs.ir_prob_scaling_factor * s_val)

    def get_immune_removal_prob(self):
        """
        Returns probability of immune cell removal due to local and global recruitment

        Probability is only non-zero if the state variable *S* is negative, in which case
        the probability is an error function of *S*

        :return: probability of immune cell removal due to local and global recruitment
        """
        s_val = self.get_state_variable_val()
        if s_val > 0:
            return 0.0
        else:
            return math.erf(- ViralInfectionVTMModelInputs.ir_prob_scaling_factor * s_val)

    def increment_total_cytokine_count(self, _inc_amount):
        """
        Increment total count of cytokine

        :param _inc_amount: total amount increment
        :return: None
        """
        self.__total_cytokine = max(0.0, self.__total_cytokine + _inc_amount)

    def set_cytokine_decay(self, _val: float):
        if _val < 0:
            raise ValueError('Cytokine decay rate must be non-negative')
        self.__ck_decay = _val


class oxidationAgentModelSteppable(ViralInfectionVTMSteppableBasePy):
    """
    Implements immune cell oxidizing agent cytotoxicity module

    This steppable uses a callback per registered cell type that is called on each cell of the registered type every
    time step to calculate interactions with the module oxidator field.

    Every secretion callback has the signature secr_func(_steppable, _cell, _mcs)
    - _steppable (oxidationAgentModelSteppable): self
    - _cell (cc3d.cpp.CompuCell.CellG): a cell
    - _mcs (int): the current simulation step

    Callbacks can be registered for a cell type with register_secretor_by_type, and unregistered for a cell type with
    unregister_secretor_by_type.

    By default, requires a diffusion field with module oxidator field name and CC3DML ids 'oxi_dc' and
    'oxi_decay' for the field diffusion and decay constants, respectively. These can be set with set_field_data.

    By default, module uninfected, infected and virus-releasing types are assigned the callback secr_func_epithelial
    and module immune cell type is assigned the callback secr_func_immune.
    """
    # todo: make oxidationAgentModelSteppable-related model parameters settable

    unique_key = ViralInfectionVTMLib.oxidation_steppable_key

    def __init__(self, frequency=1):
        super().__init__(frequency)

        self.runBeforeMCS = 1

        if ViralInfectionVTMModelInputs.track_model_variables:
            attribute_name = ViralInfectionVTMLib.oxi_killed_cellg_key
            field_name = attribute_name.replace(ViralInfectionVTMLib.module_prefix, '')
            self.track_cell_level_scalar_attribute(field_name=field_name, attribute_name=attribute_name)

        # Reference to SimDataSteppable
        self.simdata_steppable = None

        self._oxidator_diffusion_id = ''
        self._oxidator_decay_id = ''
        self._secretors_by_type = {}

        # Initialize default data
        self.set_dead_type_name(MainSteppables.dead_type_name)
        self.set_cytokine_field_name(cytokine_field_name)
        self.set_field_data(field_name=oxidator_field_name, diffusion='oxi_dc', decay='oxi_decay')
        self.register_secretor_by_type(MainSteppables.uninfected_type_name, self.secr_func_epithelial)
        self.register_secretor_by_type(MainSteppables.infected_type_name, self.secr_func_epithelial)
        self.register_secretor_by_type(MainSteppables.virus_releasing_type_name, self.secr_func_epithelial)
        self.register_secretor_by_type(immune_type_name, self.secr_func_immune)

    def start(self):
        self.get_xml_element(self._oxidator_diffusion_id).cdata = ViralInfectionVTMModelInputs.oxi_dc
        self.get_xml_element(self._oxidator_decay_id).cdata = ViralInfectionVTMModelInputs.oxi_decay

    def step(self, mcs):
        if self.simdata_steppable is None:
            try:
                self.simdata_steppable: SimDataSteppable = self.shared_steppable_vars[SimDataSteppable.unique_key]
            except KeyError:
                pass

        for type_name, secr_func in self._secretors_by_type.items():
            for cell in self.cell_list_by_type(getattr(self, type_name.upper())):
                secr_func(self, cell, mcs)

    def finish(self):
        # this function may be called at the end of simulation - used very infrequently though
        return

    def set_field_data(self, field_name: str = None, diffusion: str = None, decay: str = None):
        if field_name is not None:
            self.set_oxidator_field_name(field_name)
        if diffusion is not None:
            self._oxidator_diffusion_id = diffusion
        if decay is not None:
            self._oxidator_decay_id = decay

    def register_secretor_by_type(self, _type_name: str, _secr_func):
        """
        Register a secretion callback by cell type

        :param _type_name: name of cell type
        :param _secr_func: secretion callback
        :return: None
        """
        self._secretors_by_type[_type_name] = _secr_func

    def unregister_secretor_by_type(self, _type_name: str):
        """
        Unregister a secretion callback by cell type

        :param _type_name: name of cell type
        :return: callback of registered cell type
        """
        return self._secretors_by_type.pop(_type_name)

    @staticmethod
    def secr_func_epithelial(self, _cell, _mcs):
        seen_field = self.total_seen_field(self.oxidator_field, _cell)
        if seen_field >= ViralInfectionVTMModelInputs.oxi_death_thr:
            self.kill_cell(cell=_cell)
            _cell.dict[ViralInfectionVTMLib.oxi_killed_cellg_key] = True
            if self.simdata_steppable is not None:
                self.simdata_steppable.track_death_oxi_field()

    @staticmethod
    def secr_func_immune(self, _cell, _mcs):
        oxi_secretor = self.get_field_secretor(self._oxidator_field_name)
        if _cell.dict[ViralInfectionVTMLib.activated_cellg_key]:
            seen_field = self.total_seen_field(self.cytokine_field, _cell)
            if seen_field > ViralInfectionVTMModelInputs.oxi_sec_thr:
                oxi_secretor.secreteInsideCellTotalCount(_cell,
                                                         ViralInfectionVTMModelInputs.max_oxi_secrete / _cell.volume)
