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

from . import IFNInputs

# Module specific references
IFN_model_name = 'IFN_model'
IFN_field_name = 'IFNe'
viral_replication_model_name = 'viral_replication_model'
initial_amount_virus = 6.9e-8
min_to_mcs = s_to_mcs / 60.0
hours_to_mcs = min_to_mcs / 60.0
days_to_mcs = hours_to_mcs / 24.0


def IFN_model_string():
    """
    Antimony model string generator for IFN intracellular signaling model adopted from "Multiscale Model
    of RNA Virus Replication and Interferon Responses Reveals Factors Controlling Plaque Growth Dynamics"
    Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego, James A. Glazier and Jason E. Shoemaker
    :return: Antimony model string
    """
    model_string = f'''{IFN_model_name}()
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
        s_t = {hours_to_mcs} ;
        
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
    model_string = f'''{viral_replication_model_name}()
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
        V = 0.0      ;
        H = 1.0      ;
    
        //Inputs
        IFNe = 0.0   ;
    '''
    return model_string

ifn_model_vars = ["IFN", "STATP", "IRF7", "IRF7P"]
viral_replication_model_vars = ["H", "V"]

#TODO: Where is MainSteppables?
class IFNSteppable(ViralInfectionVTMSteppableBasePy)
    """
    Implements IFN model
    
    """
    def __init__(self, frequency=1):
        super().__init__(frequency)

        # Initialize default data
        self.register_type(MainSteppables.uninfected_type_name)
        self.register_type(MainSteppables.infected_type_name)
        self.register_type(MainSteppables.virus_releasing_type_name)

    def start(self):
        # Load IFN sbml model
        self.set_sbml_global_options(self.sbml_options)
        self.register_ode_model(model_name=IFN_model_name,
                                #TODO: what is the epithelial_model_fcn_generator()
                                model_fcn=self.epithelial_model_fcn_generator(),
                                cell_types=self._registered_types,
                                #TODO: this was a variable, but given that the model is rescaled
                                # inside the Antimony string, is reduntant
                                step_size=1.0)

        # Load viral replication sbml model
        self.register_ode_model(model_name=viral_replication_model_name,
                                #TODO: what is the epithelial_model_fcn_generator()
                                model_fcn=self.epithelial_model_fcn_generator(),
                                cell_types=self._registered_types,
                                #TODO: this was a variable, but given that the model is rescaled
                                # inside the Antimony string, is reduntant
                                step_size=1.0)

        # Set IFNe diffusion parameters
        self.get_xml_element('IFNe_dc').cdata = IFNInputs.IFNe_diffusion_coefficient * min_to_mcs
        self.get_xml_element('IFNe_decay').cdata = IFNInputs.t2 * hours_to_mcs

        # Set Virus diffusion parameters
        self.get_xml_element('virus_dc').cdata = IFNInputs.virus_diffusion_coefficient * min_to_mcs
        self.get_xml_element('virus_decay').cdata = IFNInputs.c * days_to_mcs

        # Set secretors
        self.secretorIFN = self.get_field_secretor(field_name=self.IFN_field_name)
        self.secretorV = self.get_field_secretor(field_name=self.Virus_field_name)

    def step(self,mcs):
        #TODO: Cut the repetition of the references to SBMLs

        # Production of extracellular virus
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            virus_cell_sbml = getattr(cell.sbml, viral_replication_model_name)
            k73 = virus_cell_sbml['k73']
            internalVirus = virus_cell_sbml['V']
            p = k73 * internalVirus * 1094460.28
            self.secretorV.secreteInsideCellTotalCount(cell, p / cell.volume)

        # Production of IFNe
        for cell in self.cell_list_by_type(*self._registered_types):
            IFN_cell_sbml = getattr(cell.sbml, IFN_model_name)
            k21 = IFN_cell_sbml['k21']
            intracellularIFN = IFN_cell_sbml['IFN']
            p = k21 * intracellularIFN
            self.secretorIFN.secreteInsideCellTotalCount(cell, p / cell.volume)

        # Viral cell death
        #TODO: How to substitute the cell death in the original model by this
        for cell in self.cell_list_by_type(self.VIRUSRELEASING):
            virus_cell_sbml = getattr(cell.sbml, viral_replication_model_name)
            k61 = virus_cell_sbml['k61']
            H = virus_cell_sbml['H']
            V = virus_cell_sbml['V']
            viral_death_rate = k61 * V * (1 - H)
            pr = nCoVUtils.ul_rate_to_prob(viral_death_rate)
            if random.random() <= pr:
                self.set_cell_type(cell, self.DYING)

        # Connect viral replication model and IFN model and read IFNe field
        for cell in self.cell_list_by_type(*self._registered_types):
            IFN_cell_sbml = getattr(cell.sbml, IFN_model_name)
            virus_cell_sbml = getattr(cell.sbml, viral_replication_model_name)
            IFN_cell_sbml['H'] = virus_cell_sbml['H']
            IFN_cell_sbml['V'] = virus_cell_sbml['V']
            IFN_cell_sbml['IFNe'] = self.secretorIFN.amountSeenByCell(cell)
            virus_cell_sbml['IFNe'] = self.secretorIFN.amountSeenByCell(cell)

            # Step the models for this cell
            ViralInfectionVTMLib.step_sbml_model_cell(cell=cell, name=IFN_model_name)
            ViralInfectionVTMLib.step_sbml_model_cell(cell=cell, name=viral_replication_model_name)

#TODO: Make sure that the SteppableBasePy is correct
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
            colors = ['magenta','blue', 'green', 'purple']
            for i in range(len(ifn_model_vars)):
                attr_name = 'ifn_data_win' + str(i)
                new_window = self.add_new_plot_window(title=ifn_model_vars[i],
                                                              x_axis_title='Time (hrs)',
                                                              y_axis_title='Variables', x_scale_type='linear',
                                                              y_scale_type='linear',
                                                              grid=False,
                                                              config_options={'legend': True})
                new_window.add_plot(ifn_model_vars[i], style='Dots', color=colors[i], size=5)
                setattr(self, attr_name, new_window)

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

        #Plotting and writing average values of IFN model variables
        #TODO: Instead of checking if the cell has the SBML model,
        # can I iterate over the list of types registered in the previous steppable?
        if plot_ifn_data or write_ifn_data:
            for i in range(len(ifn_model_vars)):
                L = 0.0
                total_var = 0.0
                for cell in self.cell_list_by_type(*self._registered_types):
                    if hasattr(cell.sbml, IFN_model_name):
                        cell_sbml = getattr(cell.sbml, IFN_model_name)
                        L += 1.0
                        total_var += cell_sbml[ifn_model_vars[i]]
                if plot_ifn_data:
                    window_name = getattr(self, 'ifn_data_win' + str(i))
                    window_name.add_data_point(ifn_model_vars[i], mcs * hours_to_mcs, total_var / L)
                if write_ifn_data:
                    self.ifn_data[mcs] = mcs * s_to_mcs / 60.0 / 60.0
                    self.ifn_data[mcs].append(total_var / L)

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

#TODO: ADD VIRUS REPLICATION PLOT AND WRITE, CHANGE PLOT WIN NUMBERS TO NAMES