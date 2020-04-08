# This is a general library of CompuCell3D steppable classes for the shared coronavirus modeling and simulation project
# hosted by the Biocomplexity Institute at Indiana University

from cc3d.core.PySteppables import *


class nCoVSteppableBase(SteppableBasePy):

    def __init__(self, frequency=1):
        SteppableBasePy.__init__(self, frequency)

    def start(self):
        pass

    def step(self, mcs):
        pass

    def finish(self):
        pass
