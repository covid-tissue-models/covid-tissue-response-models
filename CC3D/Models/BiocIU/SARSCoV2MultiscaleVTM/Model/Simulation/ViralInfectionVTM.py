"""
To cite this framework please use the following:

T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
Jessica Crawshaw, James M. Osborne, Ellen M. Quardokus, Richard K. Plemper, James A. Glazier,
"A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and
immune response in epithelial tissues and its application to drug therapy timing and effectiveness",
PLoS Comput Biol 16(12): e1008451. https://doi.org/10.1371/journal.pcbi.1008451
"""

import os
from ViralInfectionVTMModelInputs import __file__ as main_step_file
sys.path.append(os.path.dirname(os.path.dirname(main_step_file)))
os.environ["ViralInfectionVTM"] = os.path.dirname(os.path.dirname(main_step_file))

from cc3d import CompuCellSetup

# All imports, manipulations and registrations should occur after this line

from ViralInfectionVTMSteppables import CellInitializerSteppable
steppable = CellInitializerSteppable(frequency=1)
steppable.voxel_length = 4.0
steppable.step_period = 5.0 * 60
CompuCellSetup.register_steppable(steppable=steppable)

from ViralInfectionVTMSteppables import VirusFieldInitializerSteppable

CompuCellSetup.register_steppable(steppable=VirusFieldInitializerSteppable(frequency=1))

from Models.RecoveryNeighbor.RecoverySteppables import NeighborRecoverySteppable

CompuCellSetup.register_steppable(steppable=NeighborRecoverySteppable(frequency=1))

from ViralInfectionVTMSteppables import ViralDeathSteppable

CompuCellSetup.register_steppable(steppable=ViralDeathSteppable(frequency=1))

from ViralInfectionVTMSteppables import EclipsePhaseSteppable

CompuCellSetup.register_steppable(steppable=EclipsePhaseSteppable(frequency=1))

from ViralInfectionVTMSteppables import ViralInternalizationSteppable

CompuCellSetup.register_steppable(steppable=ViralInternalizationSteppable(frequency=1))

from ViralInfectionVTMSteppables import ViralReleaseSteppable

CompuCellSetup.register_steppable(steppable=ViralReleaseSteppable(frequency=1))

from ViralInfectionVTMSteppables import SimDataSteppable

CompuCellSetup.register_steppable(steppable=SimDataSteppable(frequency=1))

CompuCellSetup.run()
