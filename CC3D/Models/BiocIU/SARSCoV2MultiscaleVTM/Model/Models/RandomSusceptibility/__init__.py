"""
Randomly distributed virus susceptibility
=========================================

Written by T.J. Sego, Ph.D.

A fraction of uninfected cells are randomly selected and made unsusceptible to viral internalization
Model parameters are specified in SusceptibilityModelInputs.py

Maintainer(s)
=============
- T.J. Sego, Ph.D., Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

Contents
========
- SusceptibilityModelInputs.py: Default model parameters

- SusceptibilitySteppables.py: Module steppables

- example/

    - ViralInfectionVTM.py: a simulation script for running the entire module in its default configuration

    - ViralInfectionVTM.xml: a cc3dml script for running the entire module in its default configuration

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
