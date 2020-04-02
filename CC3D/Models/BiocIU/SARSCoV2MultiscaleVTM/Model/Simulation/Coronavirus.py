###############################################################################################################
# To cite this model please use the following:
#
# Josua Aponte-Serrano, T.J. Sego, James A. Glazier,
# "Model of Viral Tissue Infection"
# https://github.com/covid-tissue-models/covid-tissue-response-models/tree/master/tellurium/simple_virus_model
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

from CoronavirusSteppables import RecoverySteppable
CompuCellSetup.register_steppable(steppable=RecoverySteppable(frequency=1))

from CoronavirusSteppables import ImmuneCellSeedingSteppable
CompuCellSetup.register_steppable(steppable=ImmuneCellSeedingSteppable(frequency=1))


        
from CoronavirusSteppables import CytokineProductionAbsorptionSteppable
CompuCellSetup.register_steppable(steppable=CytokineProductionAbsorptionSteppable(frequency=1))

CompuCellSetup.run()
