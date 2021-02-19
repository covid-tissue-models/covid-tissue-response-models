"""
Data defined here corresponds to the module SimDataSteppable and framework BatchRun module
"""

from . import module_prefix

export_data_desc = {module_prefix + 'ir_data': ['ImmuneResp'],
                    module_prefix + 'med_diff_data': ['MedViral',
                                                      'MedCyt',
                                                      'MedOxi'],
                    module_prefix + 'pop_data': ['Uninfected',
                                                 'Infected',
                                                 'InfectedSecreting',
                                                 'Dying',
                                                 'ImmuneCell',
                                                 'ImmuneCellActivated'],
                    module_prefix + 'spat_data': ['DeathComp',
                                                  'InfectDist'],
                    module_prefix + 'death_data': ['Viral',
                                                   'OxiField',
                                                   'Contact',
                                                   'Bystander']}

y_label_str = {module_prefix + 'ir_data': {'ImmuneResp': 'Immune response state variable'},
               module_prefix + 'med_diff_data': {'MedViral': 'Total diffusive virus',
                                                 'MedCyt': 'Total diffusive cytokine',
                                                 'MedOxi': 'Total oxidative agent'},
               module_prefix + 'pop_data': {'Uninfected': 'Number of uninfected cells',
                                            'Infected': 'Number of infected cells',
                                            'InfectedSecreting': 'Number of infected secreting cells',
                                            'Dying': 'Number of dying cells',
                                            'ImmuneCell': 'Number of immune cells',
                                            'ImmuneCellActivated': 'Number of activated immune cells'},
               module_prefix + 'spat_data': {'DeathComp': 'Cell death compactness (ul)',
                                             'InfectDist': 'Infection distance (px)'},
               module_prefix + 'death_data': {'Viral': 'Number of virally-induced apoptosis deaths',
                                              'OxiField': 'Number of oxidative deaths',
                                              'Contact': 'Number of cytotoxic kill deaths',
                                              'Bystander': 'Number of bystander effect deaths'}
               }

fig_save_names = {module_prefix + 'ir_data': {'ImmuneResp': module_prefix + 'metric_immune_response_svar'},
                  module_prefix + 'med_diff_data': {'MedViral': module_prefix + 'metric_diffusive_virus',
                                                    'MedCyt': module_prefix + 'metric_diffusive_cytokine',
                                                    'MedOxi': module_prefix + 'metric_diffusive_oxidator'},
                  module_prefix + 'pop_data': {'Uninfected': module_prefix + 'metric_num_uninfected',
                                               'Infected': module_prefix + 'metric_num_infected',
                                               'InfectedSecreting': module_prefix + 'metric_num_infectedSecreting',
                                               'Dying': module_prefix + 'metric_num_dying',
                                               'ImmuneCell': module_prefix + 'metric_num_immune',
                                               'ImmuneCellActivated': module_prefix + 'metric_num_immuneActivated'},
                  module_prefix + 'spat_data': {'DeathComp': module_prefix + 'metric_death_compact',
                                                'InfectDist': module_prefix + 'metric_infection_distance'},
                  module_prefix + 'death_data': {'Viral': module_prefix + 'metric_death_viral',
                                                 'OxiField': module_prefix + 'metric_death_oxi',
                                                 'Contact': module_prefix + 'metric_death_contact',
                                                 'Bystander': module_prefix + 'metric_death_bystander'}
                  }
