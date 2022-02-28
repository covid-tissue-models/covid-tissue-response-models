"""
Written by JFG


"""

from numpy import arange, around

diffusion_investigation_only_change = {'media_selection': list(range(4))}

diffusion_investigation = {'media_selection': list(range(5))}

lytic_non_lytic_inv = {'lytic_prob_mult': list(around(arange(0, 1.1, 0.1),1))}

