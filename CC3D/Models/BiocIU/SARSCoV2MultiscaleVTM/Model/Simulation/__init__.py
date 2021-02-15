"""
Simple population model of infection
====================================

Written by T.J. Sego, Ph.D.

Module steppables define a spatiotemporal, multicellular implementation of the following ODE model of susceptible `S`,
exposed `E`, infectious `I`, and dead `D` cells, and viral load `V`,

dS/dt = -a * S * V

dE/dt = a * S * V - b * E

dI/dt = b * E - c * I

dD/dt = c * I

dV / dt = f * I - g * V

Default parameters are specified in ViralInfectionVTMModelInputs.py.

To cite this framework please use the following:

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
- ViralInfectionVTM.py: the cc3d execution script of this framework

- ViralInfectionVTMModelInputs.py: default parameters of this module

- ViralInfectionVTMSteppables.py: steppables

- ViralInfectionVTM.xml: cc3dml script of this framework

Change log
==========

0.0.0
-----
- Implementation as described at https://doi.org/10.1371/journal.pcbi.1008451.

"""
version_major = 0
version_minor = 0
version_build = 0
version_str = f"{version_major}.{version_minor}.{version_build}"
