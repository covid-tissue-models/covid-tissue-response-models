"""
Info for documentation of parameter values goes here, in *__param_desc__*

The parameter name and assigned description can be quickly exported to csv along with specified values using
*export_parameters*, in export_parameters.py

To use, assign the name of the code variable as a key, and a string description as its value
e.g.,
    __param_desc__['my_var'] = 'My variable'

    my_var = 1

The generated csv will read, "my_var, 1, My variable'
"""
__param_desc__ = {}

# Data control options
__param_desc__['plot_pop_data_freq'] = 'Plot population data frequency'
plot_pop_data_freq = 10  # Plot population data frequency (disable with 0)
__param_desc__['write_pop_data_freq'] = 'Write population data to simulation directory frequency'
write_pop_data_freq = 0  # Write population data to simulation directory frequency (disable with 0)
__param_desc__['plot_med_diff_data_freq'] = 'Plot total diffusive field amount frequency'
plot_med_diff_data_freq = 0  # Plot total diffusive field amount frequency (disable with 0)
__param_desc__['write_med_diff_data_freq'] = 'Write total diffusive field amount frequency'
write_med_diff_data_freq = 0  # Write total diffusive field amount frequency (disable with 0)

# Conversion Factors
__param_desc__['s_to_mcs'] = 'Simulation step'
s_to_mcs = 10 * 60  # s/mcs
__param_desc__['um_to_lat_width'] = 'Lattice width'
um_to_lat_width = 4.0  # um/lattice_length

# Experimental Parameters
__param_desc__['cell_diameter'] = 'Cell diameter'
cell_diameter = 12.0  # um

__param_desc__['secretion_rate'] = 'Secretion rate of virus'
secretion_rate = 1.0 / (24.0 * 60.0 * 60.0)  # units virus/s

__param_desc__['virus_dc'] = 'Viral diffusion coefficient'
virus_dc = 10.0 / 1000.0  # um^2/s

__param_desc__['virus_dl'] = 'Virus diffusion length'
virus_dl = cell_diameter * 3.0  # um

__param_desc__['virus_decay'] = 'Virus decay coefficient'
virus_decay = virus_dc / (virus_dl ** 2)  # 1/s

# Viral Internalization parameter
__param_desc__['internalization_rate'] = 'Internalization rate'
internalization_rate = 500.0 / (24.0 * 60.0 * 60.0)  # 1/(unit virus*s)

# Viral eclipse phase parameter
__param_desc__['eclipse_phase'] = 'Eclipse phase'
eclipse_phase = 0.25 * (24.0 * 60.0 * 60.0)  # s

# Viral death rate parameter
__param_desc__['viral_death_rate'] = 'Viral death rate parameter'
viral_death_rate = 3.0 / eclipse_phase  # 1/s

# =============================
# CompuCell Parameters
# cell
__param_desc__['volume_lm'] = 'Volume constraint LM'
volume_lm = 9
