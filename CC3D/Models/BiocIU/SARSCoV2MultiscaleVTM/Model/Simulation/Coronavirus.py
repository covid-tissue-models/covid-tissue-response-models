###############################################################################################################
# To cite this model please use the following:
#
# Josua Aponte-Serrano, T.J. Sego, Juliano F. Gianlupi, James A. Glazier,
# "Model of Viral Tissue Infection"
# https://github.com/covid-tissue-models/covid-tissue-response-models/tree/master/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM
###############################################################################################################

from cc3d import CompuCellSetup

from CoronavirusSteppables import CellsInitializerSteppable
CompuCellSetup.register_steppable(steppable=CellsInitializerSteppable(frequency=1))

from CoronavirusSteppables import Viral_ReplicationSteppable
CompuCellSetup.register_steppable(steppable=Viral_ReplicationSteppable(frequency=1))

from CoronavirusSteppables import Viral_SecretionSteppable
CompuCellSetup.register_steppable(steppable=Viral_SecretionSteppable(frequency=1))

from CoronavirusSteppables import ImmuneCellKillingSteppable
CompuCellSetup.register_steppable(steppable=ImmuneCellKillingSteppable(frequency=1))

from CoronavirusSteppables import ChemotaxisSteppable
CompuCellSetup.register_steppable(steppable=ChemotaxisSteppable(frequency=1))

from CoronavirusSteppables import ImmuneCellSeedingSteppable
CompuCellSetup.register_steppable(steppable=ImmuneCellSeedingSteppable(frequency=1))

from CoronavirusSteppables import SimDataSteppable
CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))
        
from CoronavirusSteppables import CytokineProductionAbsorptionSteppable
CompuCellSetup.register_steppable(steppable=CytokineProductionAbsorptionSteppable(frequency=1))

from CoronavirusSteppables import Viral_InternalizationSteppable
CompuCellSetup.register_steppable(steppable=Viral_InternalizationSteppable(frequency=1))

from CoronavirusSteppables import ImmuneRecruitmentSteppable
CompuCellSetup.register_steppable(steppable=ImmuneRecruitmentSteppable(frequency=1))

from CoronavirusSteppables import oxidationAgentModelSteppable
CompuCellSetup.register_steppable(steppable=oxidationAgentModelSteppable(frequency=1))

CompuCellSetup.run()
