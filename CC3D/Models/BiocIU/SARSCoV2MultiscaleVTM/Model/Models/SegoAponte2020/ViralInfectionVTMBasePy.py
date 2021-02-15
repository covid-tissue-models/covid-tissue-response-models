"""
Defines steppable base class
"""

# Import project libraries
from . import ViralInfectionVTMLib
from . import ViralInfectionVTMModelInputs

# Import toolkit
from nCoVToolkit.nCoVSteppableBase import nCoVSteppableBase


class ViralInfectionVTMSteppableBasePy(nCoVSteppableBase):
    """
    Base class of all steppables defined in this module.

    Each steppable, at most, expects the following objects with corresponding settable names as attribute

    - an uninfected cell type: ``uninfected_type_name``
    - an infected cell type: ``infected_type_name``
    - a virus-releasing cell type: ``virus_releasing_type_name``
    - a dead cell type: ``dead_type_name``
    - an immune cell type: ``immune_type_name``
    - a virus field: ``virus_field_name``
    - a cytokine field: ``cytokine_field_name``
    - an oxidative agent field: ``oxidator_field_name``

    Some steppables refer to model parameters that are defined in units of seconds.
    Conversions with these parameters are made on the fly with a settable conversion parameter ``s_to_mcs``.
    """

    vr_model_name = ViralInfectionVTMLib.vr_model_name
    _viral_replication_model_string_gen = ViralInfectionVTMLib.viral_replication_model_string

    def __init__(self, frequency=1):
        nCoVSteppableBase.__init__(self, frequency)

        self._uninfected_type_name = ''
        self._infected_type_name = ''
        self._virus_releasing_type_name = ''
        self._dead_type_name = ''
        self._immune_type_name = ''

        self._virus_field_name = ''
        self._cytokine_field_name = ''
        self._oxidator_field_name = ''

        self.s_to_mcs = ViralInfectionVTMModelInputs.s_to_mcs

    def viral_replication_model_string(self, *args, **kwargs):
        """
        Antimony model string generator for viral replication model; can be set with set_viral_replication_model(), and
        subclasses can override

        :param args:
        :param kwargs:
        :return {str}: Antimony model string for viral replication model
        """
        return ViralInfectionVTMSteppableBasePy._viral_replication_model_string_gen(*args, **kwargs)

    @staticmethod
    def set_viral_replication_model(_fnc, _name: str = ViralInfectionVTMLib.vr_model_name):
        """
        Sets the viral replication model

        :param _fnc: Antimony model string generator
        :param _name: name of Antimony model
        :return: None
        """
        ViralInfectionVTMSteppableBasePy._viral_replication_model_string_gen = _fnc
        ViralInfectionVTMSteppableBasePy.vr_model_name = _name

    def load_viral_replication_model(self, *args, **kwargs):
        """
        Loads viral replication model for a cell; subclasses can override

        :param args:
        :param kwargs:
        :return: None
        """
        assert 'cell' in kwargs, 'Specify viral replication model cell with keyword cell'
        assert 'vr_step_size' in kwargs, 'Specify viral replication model step size with keyword vr_step_size'
        return self._load_viral_replication_model(**kwargs)

    def _load_viral_replication_model(self, cell, vr_step_size, unpacking_rate=0, replicating_rate=0, r_half=0,
                                      translating_rate=0, packing_rate=0, secretion_rate=0):
        """
        Loads viral replication model for a cell; initial values of state model are extract from cell.dict

        :param cell: cell for which the viral replication model is loaded
        :param vr_step_size: Antimony/SBML model step size
        :param unpacking_rate: model unpacking rate
        :param replicating_rate: model replicating rate
        :param r_half: Value of R at which the replication rate is half max
        :param translating_rate: model translating rate
        :param packing_rate: model packing rate
        :param secretion_rate: model secretion rate
        :return: None
        """
        if cell.dict[ViralInfectionVTMLib.vrl_key]:
            self.delete_sbml_from_cell(ViralInfectionVTMSteppableBasePy.vr_model_name, cell)

        # Generate Antimony model string
        model_string = self.viral_replication_model_string(
            unpacking_rate, replicating_rate, r_half, translating_rate, packing_rate, secretion_rate,
            cell.dict[ViralInfectionVTMLib.vrm_unpacking],
            cell.dict[ViralInfectionVTMLib.vrm_replicating],
            cell.dict[ViralInfectionVTMLib.vrm_packing],
            cell.dict[ViralInfectionVTMLib.vrm_assembled],
            cell.dict[ViralInfectionVTMLib.vrm_uptake])
        self.add_antimony_to_cell(model_string=model_string,
                                  model_name=ViralInfectionVTMSteppableBasePy.vr_model_name,
                                  cell=cell,
                                  step_size=vr_step_size)
        cell.dict[ViralInfectionVTMLib.vrl_key] = True
        ViralInfectionVTMLib.enable_viral_secretion(cell, cell.type == self.virus_releasing_type_id)

    def new_cell_in_time(self, cell_type, mcs=None):
        """
        Add cell and record MCS

        :param cell_type: type id of cell (e.g., for cell.type)
        :param mcs: step when cell is created; defaults from steppable mcs attribute
        :return: new cell instance
        """
        cell = self.new_cell(cell_type)
        if mcs is None:
            if self.mcs < 0:
                mcs = 0
            else:
                mcs = self.mcs

        cell.dict[ViralInfectionVTMLib.new_cell_mcs_key] = mcs
        return cell

    def new_uninfected_cell_in_time(self, mcs=None):
        """
        Add an uninfected cell with default initial configuration and record MCS

        :param mcs: step when cell is created; defaults from steppable mcs attribute
        :return: new cell instance of immune cell type
        """
        cell = self.new_cell_in_time(self.uninfected_type_id, mcs)
        cell.dict[ViralInfectionVTMLib.vrl_key] = False
        ViralInfectionVTMLib.reset_viral_replication_variables(cell=cell)
        return cell

    def new_immune_cell_in_time(self, ck_production, ck_consumption, mcs=None, activated=False):
        """
        Add an immune cell with default initial configuration and record MCS

        :param ck_production: cytokine production rate
        :param ck_consumption: cytokine consumption rate
        :param mcs: step when cell is created; defaults from steppable mcs attribute
        :param activated: flag for immune cell being naive or activated
        :return: new cell instance of immune cell type
        """
        cell = self.new_cell_in_time(self.immune_type_id, mcs)
        # cyttokine params
        cell.dict[ViralInfectionVTMLib.ck_production_cellg_key] = ck_production
        cell.dict[ViralInfectionVTMLib.ck_consumption_cellg_key] = ck_consumption
        cell.dict[ViralInfectionVTMLib.activated_cellg_key] = activated
        cell.dict[ViralInfectionVTMLib.tot_ck_upt_cellg_key] = 0
        return cell

    def total_seen_field(self, field, cell, estimate=True):
        """
        Calculates total value of field in the cell.

        :param field: the field to be looked
        :param cell: the cell
        :param estimate: when true assumes homogeneous field. false is slower
        :return: calculated total field value
        """
        if estimate:
            tot_field = field[cell.xCOM, cell.yCOM, cell.zCOM] * cell.volume
        else:
            tot_field = 0
            for ptd in self.get_cell_pixel_list(cell):
                tot_field += field[ptd.pixel.x, ptd.pixel.y, ptd.pixel.z]

        return tot_field

    def kill_cell(self, cell):
        """
        Model-specific cell death routines

        :param cell: cell to kill
        :return: None
        """
        self.set_cell_type(cell, self.dead_type_id)

        # Remove viral replication model: no model for dead cell type
        ViralInfectionVTMLib.reset_viral_replication_variables(cell=cell)
        self.remove_viral_replication_model(cell=cell)

    def remove_viral_replication_model(self, cell):
        """
        Removes viral replication model for a cell

        :param cell: cell for which to remove the viral replication model
        :return: None
        """
        if cell.dict[ViralInfectionVTMLib.vrl_key]:
            self.delete_sbml_from_cell(ViralInfectionVTMSteppableBasePy.vr_model_name, cell)
            cell.dict[ViralInfectionVTMLib.vrl_key] = False

    def set_uninfected_type_name(self, _name: str):
        """
        Set the name of the uninfected cell type according to a cc3d simulation

        :param _name: Name of the uninfected cell type according to a cc3d simulation
        :return: None
        """
        self._uninfected_type_name = _name

    @property
    def uninfected_type_name(self):
        """
        Name of the uninfected cell type according to a cc3d simulation
        """
        return self._uninfected_type_name

    @uninfected_type_name.setter
    def uninfected_type_name(self, _name: str):
        self.set_uninfected_type_name(_name)

    def set_infected_type_name(self, _name: str):
        """
        Set the name of the infected cell type according to a cc3d simulation

        :param _name: Name of the infected cell type according to a cc3d simulation
        :return: None
        """
        self._infected_type_name = _name

    @property
    def infected_type_name(self):
        """
        Name of the infected cell type according to a cc3d simulation
        """
        return self._infected_type_name

    @infected_type_name.setter
    def infected_type_name(self, _name: str):
        self.set_infected_type_name(_name)

    def set_virus_releasing_type_name(self, _name: str):
        """
        Set the name of the virus-releasing cell type according to a cc3d simulation

        :param _name: Name of the virus-releasing cell type according to a cc3d simulation
        :return: None
        """
        self._virus_releasing_type_name = _name

    @property
    def virus_releasing_type_name(self):
        """
        Name of the virus-releasing cell type according to a cc3d simulation
        """
        return self._virus_releasing_type_name

    @virus_releasing_type_name.setter
    def virus_releasing_type_name(self, _name: str):
        self.set_virus_releasing_type_name(_name)

    def set_dead_type_name(self, _name: str):
        """
        Set the name of the dead cell type according to a cc3d simulation

        :param _name: Name of the dead cell type according to a cc3d simulation
        :return: None
        """
        self._dead_type_name = _name

    @property
    def dead_type_name(self):
        """
        Name of the dead cell type according to a cc3d simulation
        """
        return self._dead_type_name

    @dead_type_name.setter
    def dead_type_name(self, _name: str):
        self.set_dead_type_name(_name)

    def set_immune_type_name(self, _name: str):
        """
        Set the name of the immune cell type according to a cc3d simulation

        :param _name: Name of the immune cell type according to a cc3d simulation
        :return: None
        """
        self._immune_type_name = _name

    @property
    def immune_type_name(self):
        """
        Name of the immune cell type according to a cc3d simulation
        """
        return self._immune_type_name

    @immune_type_name.setter
    def immune_type_name(self, _name: str):
        self.set_immune_type_name(_name)

    @property
    def uninfected_type_id(self) -> int:
        """
        Id of the uninfected cell type according to a cc3d simulation
        """
        return getattr(self, self._uninfected_type_name.upper())

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

    @property
    def immune_type_id(self) -> int:
        """
        Id of the immune cell type according to a cc3d simulation
        """
        return getattr(self, self._immune_type_name.upper())

    def set_virus_field_name(self, _name: str):
        """
        Set the virus field name

        :param _name: Virus field name
        :return: None
        """
        self._virus_field_name = _name

    @property
    def virus_field_name(self):
        """
        Virus field name
        """
        return self._virus_field_name

    @virus_field_name.setter
    def virus_field_name(self, _name: str):
        self.set_virus_field_name(_name)

    def set_cytokine_field_name(self, _name: str):
        """
        Set the cytokine field name

        :param _name: Cytokine field name
        :return: None
        """
        self._cytokine_field_name = _name

    @property
    def cytokine_field_name(self):
        """
        Cytokine field name
        """
        return self._cytokine_field_name

    @cytokine_field_name.setter
    def cytokine_field_name(self, _name: str):
        self.set_cytokine_field_name(_name)

    def set_oxidator_field_name(self, _name: str):
        """
        Set the oxidative agent field name

        :param _name: Oxidative agent field name
        :return: None
        """
        self._oxidator_field_name = _name

    @property
    def oxidator_field_name(self):
        """
        Oxidative agent field name
        """
        return self._oxidator_field_name

    @oxidator_field_name.setter
    def oxidator_field_name(self, _name: str):
        self.set_oxidator_field_name(_name)

    @property
    def virus_field(self):
        """
        Reference to the virus field
        """
        if self._virus_field_name:
            return getattr(self.field, self._virus_field_name)
        return None

    @property
    def cytokine_field(self):
        """
        Reference to the cytokine field
        """
        if self._cytokine_field_name:
            return getattr(self.field, self._cytokine_field_name)
        return None

    @property
    def oxidator_field(self):
        """
        Reference to the oxidative agent field
        """
        if self._oxidator_field_name:
            return getattr(self.field, self._oxidator_field_name)
        return None

    def set_time_step(self, _val: float):
        """
        Set the time step, in units seconds per step

        :param _val: Time step
        :return: None
        """
        if _val <= 0.0:
            raise ValueError("Step period must be positive")
        self.s_to_mcs = _val

    def get_time_step(self) -> float:
        """
        Time step, in units seconds per step
        """
        return self.s_to_mcs
