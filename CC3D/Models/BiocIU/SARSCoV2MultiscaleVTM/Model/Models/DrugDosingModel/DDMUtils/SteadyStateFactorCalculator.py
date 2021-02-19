from Models.DrugDosingModel.rate_sets import set2 as s


def get_factor(k0=s.k0, kE0=s.kE0, kp=s.kp, kpp=s.kpp, k12=s.k12, kE1=s.kE1, k23=s.k23, kE2=s.kE2, k34=s.k34, kE3=s.kE3,
               kE4=s.kE4):
    """

    :param k0:
    :param kE0:
    :param kp:
    :param kpp:
    :param k12:
    :param kE1:
    :param k23:
    :param kE2:
    :param k34:
    :param kE3:
    :param kE4:
    :return:
    """
    term1 = k34/kE4
    term2 = k23/(k34+kE3)
    term3 = k12/(k23+kE2)
    term4 = k0/(k0+k12+kE2)
    term5 = 1/(k0+kE0+kp)
    term6 = kp + (k0*k0/(k0+k12+kE1))

    factor = term1 * term2 * term3 * term4 * term5 * term6

    return factor


if __name__ == '__main__':
    print(get_factor())
