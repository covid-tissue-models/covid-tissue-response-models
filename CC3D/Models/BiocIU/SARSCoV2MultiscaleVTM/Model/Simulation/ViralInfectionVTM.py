###############################################################################################################
# To cite this model please use the following:
#
# T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
# James M. Osborne, Ellen M. Quardokus, Richard K. Plemper, James A. Glazier,
# "A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and
# immune response in epithelial tissues and its application to drug therapy timing and effectiveness",
# bioRxiv 2020.04.27.064139
###############################################################################################################

import os
from ViralInfectionVTMSteppables import __file__ as main_step_file
sys.path.append(os.path.dirname(os.path.dirname(main_step_file)))
os.environ["ViralInfectionVTM"] = os.path.dirname(os.path.dirname(main_step_file))

from cc3d import CompuCellSetup

from ViralInfectionVTMSteppables import CellsInitializerSteppable

CompuCellSetup.register_steppable(steppable=CellsInitializerSteppable(frequency=1))

from ViralInfectionVTMSteppables import ViralInternalizationSteppable

CompuCellSetup.register_steppable(steppable=ViralInternalizationSteppable(frequency=1))

from ViralInfectionVTMSteppables import ViralReplicationSteppable

CompuCellSetup.register_steppable(steppable=ViralReplicationSteppable(frequency=1))

from ViralInfectionVTMSteppables import ViralSecretionSteppable

CompuCellSetup.register_steppable(steppable=ViralSecretionSteppable(frequency=1))

from ViralInfectionVTMSteppables import ImmuneCellKillingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellKillingSteppable(frequency=1))

from ViralInfectionVTMSteppables import ChemotaxisSteppable

CompuCellSetup.register_steppable(steppable=ChemotaxisSteppable(frequency=1))

from ViralInfectionVTMSteppables import ImmuneCellSeedingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellSeedingSteppable(frequency=1))

from ViralInfectionVTMSteppables import SimDataSteppable

CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))

from ViralInfectionVTMSteppables import CytokineProductionAbsorptionSteppable

CompuCellSetup.register_steppable(steppable=CytokineProductionAbsorptionSteppable(frequency=1))

from ViralInfectionVTMSteppables import ImmuneRecruitmentSteppable

CompuCellSetup.register_steppable(steppable=ImmuneRecruitmentSteppable(frequency=1))

from ViralInfectionVTMSteppables import oxidationAgentModelSteppable

CompuCellSetup.register_steppable(steppable=oxidationAgentModelSteppable(frequency=1))

from Models.DrugDosingModel.DrugDosingModelSteppable import DrugDosingModelSteppable

CompuCellSetup.register_steppable(steppable=DrugDosingModelSteppable(frequency=1))

from Models.DrugDosingModel.DrugDosingModelSteppable import DrugDosingDataFieldsPlots
CompuCellSetup.register_steppable(steppable=DrugDosingDataFieldsPlots(frequency=1))
#
# from Models.DrugDosingModel.DrugDosingModelSteppable import ProdrugDiffusionController
# CompuCellSetup.register_steppable(steppable=ProdrugDiffusionController(frequency=1))

CompuCellSetup.run()
