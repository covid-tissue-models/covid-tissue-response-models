"""
Model of neighbor-dependent recovery
====================================

Written by T.J. Sego, Ph.D., and presented in the Interactive Two-Part Virtual Mini-workshop on
Open-Source CompuCell3D Multiscale, Virtual-Tissue Spatiotemporal Modeling and Simulations of COVID-19 Infection,
Viral Spread and Immune Response and Treatment Regimes, June 11-12 & Jun 18-19, 2020.

Dead cells "resurrect" and become uninfected according to the number of uninfected cells in its neighborhood and a
rate "recovery_rate" defined in RecoveryInputs.py.

Built from RecoverySimple module written by T.J. Sego, Ph.D.

Maintainer(s)
=============
- T.J. Sego, Ph.D., Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

Contents
========
- RecoveryInputs.py: Default model parameters

- RecoverySteppables.py: Module steppables

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
