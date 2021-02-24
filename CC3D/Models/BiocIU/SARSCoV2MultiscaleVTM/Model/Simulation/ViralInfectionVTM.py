###############################################################################################################
# To cite this model please use the following:
#
# T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
# Jessica Crawshaw, James M. Osborne, Ellen M. Quardokus, Richard K. Plemper, James A. Glazier,
# "A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and
# immune response in epithelial tissues and its application to drug therapy timing and effectiveness",
# PLoS Comput Biol 16(12): e1008451. https://doi.org/10.1371/journal.pcbi.1008451
###############################################################################################################

import os
from ViralInfectionVTMModelInputs import __file__ as main_step_file
sys.path.append(os.path.dirname(os.path.dirname(main_step_file)))
os.environ["ViralInfectionVTM"] = os.path.dirname(os.path.dirname(main_step_file))

from cc3d import CompuCellSetup

from ..Models.IFNSIgnaling.IFNSteppables import IFNVirusFieldInitializerSteppable

CompuCellSetup.register_steppable(steppable=IFNVirusFieldInitializerSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNViralDeathSteppable

CompuCellSetup.register_steppable(steppable=IFNViralDeathSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNEclipsePhaseSteppable

CompuCellSetup.register_steppable(steppable=IFNEclipsePhaseSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNViralInternalizationSteppable

CompuCellSetup.register_steppable(steppable=IFNViralInternalizationSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNViralReleaseSteppable

CompuCellSetup.register_steppable(steppable=IFNViralReleaseSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNReleaseSteppable

CompuCellSetup.register_steppable(steppable=IFNReleaseSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNFieldInitializerSteppable

CompuCellSetup.register_steppable(steppable=IFNFieldInitializerSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNSimDataSteppable

CompuCellSetup.register_steppable(steppable=IFNSimDataSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNSteppable

CompuCellSetup.register_steppable(steppable=IFNSteppable(frequency=1))

from ..Models.IFNSIgnaling.IFNSteppables import IFNCellInitializerSteppable

CompuCellSetup.register_steppable(steppable=IFNCellInitializerSteppable(frequency=1))

CompuCellSetup.run()
