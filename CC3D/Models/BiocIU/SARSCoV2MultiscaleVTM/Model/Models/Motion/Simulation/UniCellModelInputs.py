# Info for documentation of parameter values goes here, in *__param_desc__*
# The parameter name and assigned description can be quickly exported to csv along with specified values using
# *export_parameters*, in export_parameters.py
# To use, assign the name of the code variable as a key, and a string description as its value
# e.g.,     __param_desc__['my_var'] = 'My variable'
#           my_var = 1
# The generated csv will read, "my_var, 1, My variable'
__param_desc__ = {}

# Data control options
__param_desc__['track_model_variables'] = 'Enables cell-level tracking of model variables (for rendering)'
track_model_variables = False  # Set to true to enable cell-level tracking of model variables (for rendering)
__param_desc__['plot_vrm_data_freq'] = 'Plot viral replication model data frequency'
plot_vrm_data_freq = 0  # Plot viral replication model data frequency (disable with 0)
__param_desc__['write_vrm_data_freq'] = 'Write viral replication model data to simulation directory frequency'
write_vrm_data_freq = 0  # Write viral replication model data to simulation directory frequency (disable with 0)

__param_desc__["alpha"] = "persistencia"
alpha = 1

__param_desc__["beta"] = "ruido"
beta = 1

__param_desc__["numero"] = "numero"
numero = 1
