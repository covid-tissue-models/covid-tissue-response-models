"""
Defines default module model parameters
"""
__param_desc__ = {}

# Data control options
__param_desc__['plot_rec_data_freq'] = 'Plot recovery model data frequency'
plot_rec_data_freq = 0  # Plot recovery model data frequency (disable with 0)
__param_desc__['write_rec_data_freq'] = 'Write recovery model data to simulation directory frequency'
write_rec_data_freq = 0  # Write recovery model data to simulation directory frequency (disable with 0)

# Model inputs
__param_desc__['recovery_rate'] = "Recovery rate"
recovery_rate = 1.0 / (7.0 * 24.0 * 60.0 * 60.0)
