"""
Model of acute primary viral infection and immune response in epithelial tissues
================================================================================

Written by T.J. Sego, Ph.D., Josua Aponte-Serrano and Juliano Gianlupi, with assistance from Kira Breithaupt,
Samuel Heaps, Jairaj Mathur and James Glazier, Ph.D.

For details of specific models of this module, see https://doi.org/10.1371/journal.pcbi.1008451.

To cite this model please use the following:

T.J. Sego, Josua O. Aponte-Serrano, Juliano Ferrari Gianlupi, Samuel R. Heaps, Kira Breithaupt, Lutz Brusch,
Jessica Crawshaw, James M. Osborne, Ellen M. Quardokus, Richard K. Plemper, James A. Glazier,
"A modular framework for multiscale, multicellular, spatiotemporal modeling of acute primary viral infection and
immune response in epithelial tissues and its application to drug therapy timing and effectiveness",
PLoS Comput Biol 16(12): e1008451. https://doi.org/10.1371/journal.pcbi.1008451

Maintainer(s)
=============
- T.J. Sego, Ph.D., Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

- Josua Aponte-Serrano, Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

- Juliano Gianlupi, Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

Contents
========
- ViralInfectionVTMBasePy.py: defines a base class defining common data and methods to all steppables

- ViralInfectionVTMLib.py: defines basic definitions and methods

- ViralInfectionVTMModelInputs.py: defines default model parameters

- ViralInfectionVTMSteppables.py: defines all steppables implementing all model modules

- example/

    - ViralInfectionVTM.py: a simulation script for running the entire module in its default configuration

    - ViralInfectionVTM.xml: a cc3dml script for running the entire module in its default configuration

Change log
==========
1.0.0
-----
- Deployment per architecture defined in framework v1.0.0

0.0.0
-----
- Implementation as described at https://doi.org/10.1371/journal.pcbi.1008451.

"""

module_prefix = "vivtm_"

version_major = 1
version_minor = 0
version_build = 0
version_str = f"{version_major}.{version_minor}.{version_build}"
