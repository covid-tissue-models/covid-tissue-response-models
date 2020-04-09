# This is a general library for the shared coronavirus modeling and simulation project
# hosted by the Biocomplexity Institute at Indiana University


def hill_equation(val, diss_cf, hill_cf):
    """
    Hill equation
    :param val: input value
    :param diss_cf: dissociation coefficient
    :param hill_cf: Hill coefficient
    :return: Hill equation for input *val*
    """
    return 1 / (1 + (diss_cf / val) ** hill_cf)
