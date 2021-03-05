"""
Model of interferon signaling response and RNA viral replication in epithelial tissues
================================================================================

Written by Josua Aponte-Serrano and T.J. Sego, Ph.D., with assistance from Jordan Weaver, Jason Shoemaker, Ph. D.
and James Glazier, Ph.D.

To cite this model please use the following:

Multiscale Model of RNA Virus Replication and Interferon Responses Reveals Factors Controlling Plaque
Growth Dynamics" Josua O. Aponte-Serrano, Jordan J.A. Weaver, T.J. Sego, James A. Glazier and
Jason E. Shoemaker

Maintainer(s)
=============
- Josua Aponte-Serrano, Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

- T.J. Sego, Ph.D., Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

Contents
========

- IFNInputs.py: defines default model parameters

- IFNSteppables.py: defines steppables implementing model-specific modules

- example/

    - ViralInfectionVTM.py: a simulation script for running the entire module in its default configuration

    - ViralInfectionVTM.xml: a cc3dml script for running the entire module in its default configuration

- validation/

    - ifn_data.dat: ifn model data from representative cc3d simulation for validation purposes

    - ode_data.txt: ODE model data for validation purposes

    - plot_validation.py: python script to plot cc3d simulation data and ODE data for validation purposes.

    - validation.pdf: plots of ifn model data, reproduces figure 2 of the reference publication.

Change log
==========
1.0.0
-----
- Implementation as described in reference publication

"""

module_prefix = "ifnsig_"

version_major = 1
version_minor = 0
version_build = 0
version_str = f"{version_major}.{version_minor}.{version_build}"
