from Models.DrugDosingModel.rate_sets import set1
from Models.DrugDosingModel.rate_sets import set2
from Models.DrugDosingModel.rate_sets import set3
from Models.DrugDosingModel.rate_sets import set4
from Models.DrugDosingModel.rate_sets import set5
from Models.DrugDosingModel.rate_sets import set6
from Models.DrugDosingModel.rate_sets import set7
from Models.DrugDosingModel.rate_sets import set8
from Models.DrugDosingModel.rate_sets import set9
from Models.DrugDosingModel.rate_sets import set10


def import_sets_as_dict():
    rate_sets = {'set1': set1, 'set2': set2, 'set3': set3, 'set4': set4, 'set5': set5, 'set6': set6, 'set7': set7,
                 'set8': set8, 'set9': set9, 'set10': set10}

    return rate_sets
