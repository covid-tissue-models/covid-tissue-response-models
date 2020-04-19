# This is a general library for the shared coronavirus modeling and simulation project
# hosted by the Biocomplexity Institute at Indiana University


def export_parameters(param_module, export_file):
    """
    Exports parameters defined in a module to csv
    :param param_module: module containing parameters
    :param export_file: location of file to write csv
    :return: None
    """
    import csv

    params_list = [x for x in dir(param_module) if not (x.startswith('__') and x.endswith('__'))]
    params_dict = dict()
    for k in params_list:
        params_dict[k] = getattr(param_module, k)

    with open(export_file, 'w', newline='') as fout:
        csv_data_writer = csv.writer(fout, delimiter=',')
        [csv_data_writer.writerow([k, v]) for k, v in params_dict.items()]


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
