"""
Written by JFG

This file contains all the parameters used for the investigation of antiviral treatment of
infection by SARS-2. It is part of publication TBD

Each dictionary is a series of multipliers that get applied to the parametrs set in
cellular-model/Models/DrugDosingModel/DrugDosingInputs.py. Because of nuances of the code there are some
parameters that need to be changed directly in DrugDosingInputs.py. The most important one
being the flag for doing intercell variation of metabolization rates. To simulate that situation you can use any
of the dictionaries here and set ```intercell_var``` to True in that file.

How the scripts vary the values of the dictionaries is detailed in cellular-model/batch_run.py

"""

variation_in_sego_model = {'first_dose': [0],
                           'daily_dose': [0],
                           'dose_interval': [1],
                           'ic50_multiplier': [1],
                           'kon': [1]}

macro_variation = {'first_dose': [0],
                   'dose_interval': [8/24, 12/24, 1, 2, 3],
                   'ic50_multiplier': [0.01, 0.05, 0.1, 0.5, 1, 5, 10]}

treatment_starts_0 = {'first_dose': [0],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1]}
treatment_starts_1 = {'first_dose': [1],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1]}
treatment_starts_3 = {'first_dose': [3],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1]}

treatment_starts_0_halved_half_life = {'first_dose': [0],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1/2]}
treatment_starts_1_halved_half_life = {'first_dose': [1],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1/2]}
treatment_starts_3_halved_half_life = {'first_dose': [1],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1/2]}

treatment_starts_0_quartered_half_life = {'first_dose': [0],  #
                        'dose_interval': [1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6],
                        'ic50_multiplier': [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1],
                        't_half_mult': [1/4]}
