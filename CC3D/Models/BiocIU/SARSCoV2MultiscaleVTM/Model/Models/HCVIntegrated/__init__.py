"""
Integrated model of hepatitis c virus
=====================================

Written by T.J. Sego, Ph.D.

Genomic replication in the module SegoAponte2020 is replaced by a model of hepatitis c virus replication from

  Dahari, Harel, et al. "Mathematical modeling of subgenomic hepatitis C virus replication in Huh-7 cells."
  Journal of virology 81.2 (2007): 750-760.

Maintainer(s)
=============
- T.J. Sego, Ph.D., Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

Contents
========
- HCVInputs.py: Default model parameters

- HCVSteppables.py: Module steppables

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
version_major = 1
version_minor = 0
version_build = 0
version_str = f"{version_major}.{version_minor}.{version_build}"
