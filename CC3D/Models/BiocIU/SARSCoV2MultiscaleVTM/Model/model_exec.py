import argparse
import json
import os
import sys

from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller


STATUS_FILE = 'running.status'


class ArgParser:
    def __init__(self, argsv):
        self.parser = argparse.ArgumentParser()

        self.parser.add_argument('--input', '-i', type=str, nargs=1, dest='input',
                                 help='path to input json', required=True)
        self.parser.add_argument('--status', '-s', type=str, nargs=1, dest='status',
                                 help='optional path to status file', default=None)
        self.parser.add_argument('-g', action='store_true', dest='gen_status')
        self.parsed_args = self.parser.parse_args(argsv)

    @property
    def input_file(self) -> str:
        return self.parsed_args.input[0]

    @property
    def status_file(self):
        if self.parsed_args.status is None:
            return None
        else:
            return self.parsed_args.status[0]

    @property
    def generate_status(self) -> bool:
        return self.parsed_args.gen_status


def parse_args(argsv):
    return ArgParser(argsv)


def main(input_dict: dict, status_filename: str = None, generate_status: bool = False):
    if status_filename is not None:
        generate_status = False

    if generate_status:
        status_filename = os.path.join(input_dict['output_dir'], STATUS_FILE)
        with open(status_filename, "w"):
            pass
    simulation_fname = r"/N/u/jferrari/Carbonate/covid_models/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Model/Models/Motion/UniCell_Explicit_Force_2D_with_beta.cc3d"

    cc3d_caller = CC3DCaller(cc3d_sim_fname=simulation_fname,
                             output_frequency=input_dict['output_frequency'],
                             screenshot_output_frequency=input_dict['screenshot_output_frequency'],
                             output_dir=input_dict['output_dir'],
                             sim_input=input_dict['sim_input'])
    ret_val = cc3d_caller.run()

    if status_filename is not None and os.path.isfile(status_filename):
        os.remove(status_filename)

    return ret_val


if __name__ == '__main__':
    parser = parse_args(sys.argv[1:])

    with open(parser.input_file, 'r') as sf:
        _input_dict = json.load(sf)

    main(_input_dict, status_filename=parser.status_file, generate_status=parser.generate_status)
