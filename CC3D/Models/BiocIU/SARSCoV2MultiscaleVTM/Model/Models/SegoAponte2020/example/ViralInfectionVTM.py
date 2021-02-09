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
from ViralInfectionVTMModelInputs import __file__ as f
sys.path.append(os.path.dirname(os.path.dirname(f)))
os.environ["ViralInfectionVTM"] = os.path.dirname(os.path.dirname(f))

from cc3d import CompuCellSetup

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ViralReplicationSteppable

CompuCellSetup.register_steppable(steppable=ViralReplicationSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import CellsInitializerSteppable

CompuCellSetup.register_steppable(steppable=CellsInitializerSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ViralInternalizationSteppable

CompuCellSetup.register_steppable(steppable=ViralInternalizationSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ViralSecretionSteppable

CompuCellSetup.register_steppable(steppable=ViralSecretionSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ImmuneCellKillingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellKillingSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ChemotaxisSteppable

CompuCellSetup.register_steppable(steppable=ChemotaxisSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ImmuneCellSeedingSteppable

CompuCellSetup.register_steppable(steppable=ImmuneCellSeedingSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import SimDataSteppable

CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import CytokineProductionAbsorptionSteppable

CompuCellSetup.register_steppable(steppable=CytokineProductionAbsorptionSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import ImmuneRecruitmentSteppable

CompuCellSetup.register_steppable(steppable=ImmuneRecruitmentSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import oxidationAgentModelSteppable

CompuCellSetup.register_steppable(steppable=oxidationAgentModelSteppable(frequency=1))

from Models.SegoAponte2020.ViralInfectionVTMSteppables import VirusFieldInitializerSteppable

CompuCellSetup.register_steppable(steppable=VirusFieldInitializerSteppable(frequency=1))

CompuCellSetup.run()
