# This is a general library for the shared coronavirus modeling and simulation project
# hosted by the Biocomplexity Institute at Indiana University


def export_parameters(param_module, export_file):
    """
    Exports parameters defined in a module or dictionary to csv
    :param param_module: module or dictionary containing parameters
    :param export_file: location of file to write csv
    :return: None
    """
    import csv

    param_desc = None
    desc_key = '__param_desc__'

    if isinstance(param_module, dict):
        params_dict = param_module
        if desc_key in params_dict.keys():
            param_desc = params_dict.pop(desc_key)
    else:
        params_dict = dict()
        params_list = [x for x in dir(param_module) if not (x.startswith('__') and x.endswith('__'))]
        for k in params_list:
            params_dict[k] = getattr(param_module, k)
        if desc_key in dir(param_module):
            param_desc = getattr(param_module, desc_key)

    with open(export_file, 'w', newline='') as fout:
        csv_data_writer = csv.writer(fout, delimiter=',')
        for k, v in params_dict.items():
            if param_desc is not None:
                this_desc = ' '
                if k in param_desc.keys():
                    this_desc = param_desc[k]
                row_to_write = [k, v, this_desc]
            else:
                row_to_write = [k, v]
            csv_data_writer.writerow(row_to_write)

def hill_equation(val, diss_cf, hill_cf):
    """
    Hill equation
    :param val: input value
    :param diss_cf: dissociation coefficient
    :param hill_cf: Hill coefficient
    :return: Hill equation for input *val*
    """
    if val == 0:
        return 0
    else:
        return 1 / (1 + (diss_cf / val) ** hill_cf)
