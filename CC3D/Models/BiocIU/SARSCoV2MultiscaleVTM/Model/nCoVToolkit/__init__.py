"""
Framework Toolkit
=================

This module defines the basic functionality and utilities of the framework.

Written by T.J. Sego, Ph.D.

Maintainer(s)
=============
- T.J. Sego, Ph.D., Biocomplexity Institute, Indiana University, Bloomington, IN, U.S.A.

Contents
========
- nCoVSteppableBase.py: defines the CompuCell3D base steppable class for the framework

- nCoVUtils.py: defines basic utility functions of the framework

Change log
==========
1.0.0
-----
- Deployment per architecture defined in framework v1.0.0

0.0.0
-----
- Implementation as deployed to perform simulations shown at https://doi.org/10.1371/journal.pcbi.1008451.

"""

__all__ = ["nCoVSteppableBase",
           "nCoVUtils"]

version_major = 1
version_minor = 0
version_build = 0
version_str = f"{version_major}.{version_minor}.{version_build}"
