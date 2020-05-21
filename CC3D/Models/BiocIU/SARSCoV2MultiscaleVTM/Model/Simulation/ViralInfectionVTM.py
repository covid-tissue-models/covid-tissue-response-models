###############################################################################################################
# To cite this model please use the following:
#
# T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
# James M. Osborne, Ellen M. Quardokus, James A. Glazier,
# "A modular framework for multiscale multicellular spatial modeling of viral infection, immune response and drug
# therapy timing and efficacy in epithelial tissues",
# bioRxiv 2020.04.27.064139
###############################################################################################################

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

CompuCellSetup.run()
