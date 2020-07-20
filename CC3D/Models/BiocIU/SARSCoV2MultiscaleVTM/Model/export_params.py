from os.path import abspath, dirname, join

export_file = "ViralInfectionVTMModelParams.csv"


def export_parameters(file=export_file):
    from nCoVToolkit import nCoVUtils
    from Simulation import ViralInfectionVTMModelInputs

    export_file_abs = join(abspath(dirname(__file__)), file)
    nCoVUtils.export_parameters(ViralInfectionVTMModelInputs, export_file_abs)


if __name__ == "__main__":
    export_parameters()
