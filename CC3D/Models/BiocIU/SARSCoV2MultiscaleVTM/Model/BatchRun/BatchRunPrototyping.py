import json
import os
import random
import sys

from . import BatchRunLib
from .CallableCoV2VTM import CoV2VTMSimRunAsync, CallableCoV2VTMScheduler, simulation_fname

# Set this to False to not keep output of Carbonate
keep_output = True

# Status control is
#   Dispatcher checks for the existence of a file "running.status" in the output directory of a run before execution
#       If "running.status" is found in the output directory, it is interpreted as a previous run that failed
#   Dispatcher issues creation a file "running.status" in the output directory of a run before execution
#   The Python execution script will consume running.status after execution if told to do so

# Dispatch workflow is
#   Create the output directory according to the standard data structure
#   Generate launch .sh script
#   Generate the .job script that points to generated launch script and then call it through a generic shell script

cc3d_input_dict = {
    'simulation_fname': None,
    'output_frequency': None,
    'screenshot_output_frequency': None,
    'output_dir': None,
    'sim_input': None
}

# Per general compute nodes
_ciu_lims = {'nn': 72,
             'ppn': 24,
             'vmem': 503,
             'wh': None}

_default_config = {'ec': None,
                   'ee': None,
                   'jn': "MYJOB",
                   'ko': False,
                   'nn': 1,
                   'ppn': 1,
                   'shell_script': None,
                   'vmem': None,
                   'wh': None,
                   'wm': None}


def carbonate_config_template():
    return _default_config.copy()


class CallableCC3DCarbonateDispatcher:
    def __init__(self, model_config: dict, carb_config: dict = None):
        if not os.path.isdir(self.__work_dir):
            os.mkdir(self.__work_dir)

        self.rand_id = random.randint(0, sys.maxsize)

        # CC3D inputs
        self.__cc3d_input_dict = cc3d_input_dict.copy()
        for k, v in model_config.items():
            if k in list(cc3d_input_dict.keys()):
                self.__cc3d_input_dict[k] = v
        for k in ['simulation_fname', 'output_frequency', 'screenshot_output_frequency', 'output_dir']:
            assert self.__cc3d_input_dict[k] is not None

        if not os.path.isdir(self.__cc3d_input_dict['output_dir']):
            os.mkdir(self.__cc3d_input_dict['output_dir'])

        # Carbonate job description
        self.__config_dict = _default_config.copy()
        if carb_config is not None:
            for k, v in carb_config.items():
                if k in list(self.__config_dict.keys()):
                    self.__config_dict[k] = v

    @property
    def __work_dir(self):
        return os.path.join(os.path.dirname(__file__), "tmp")

    @property
    def json_filename(self):
        return os.path.join(self.__work_dir, f"model_input_{self.rand_id}.json")

    @property
    def shell_script_filename(self):
        return os.path.join(self.__work_dir, f"job_{self.rand_id}.sh")

    @property
    def job_script_filename(self):
        return os.path.join(self.__work_dir, f"job_{self.rand_id}.job")

    @staticmethod
    def _write_linux(filename: str, filecontents: str):
        with open(filename, "w", newline="\n") as f:
            f.write(filecontents)
        return

    def generate_input_json(self):
        """
        Writes a json with input variables and values, the path of which is passed to the Python execution script
        :return: {str} Path to the input json
        """
        print(f'Generating input json: {self.json_filename}')
        with open(self.json_filename, "w") as f:
            json.dump(self.__cc3d_input_dict, f)
        return

    def generate_shell_script(self):
        """
        Writes a shell script that, when called, will execute the Python execution script
        :return: {str} Path to the execution shell script
        """
        print(f'Generating execution shell: {self.shell_script_filename}')
        run_carbonate = f"{os.path.dirname(__file__)}/run_carbonate.sh"
        run_model = f"{os.path.dirname(os.path.dirname(__file__))}/model_exec.py"
        model_inputs = f"{self.json_filename}"
        o = "#!/bin/bash\n"
        o += f"sh {run_carbonate} -r {run_model} -i {model_inputs} -g\n"
        self._write_linux(self.shell_script_filename, o)

    def generate_job_script(self):
        """
        Writes a string for a .job script that, when called, will execute the shell script
        :return: {str} Path to the .job script
        """
        print(f'Generating job script: {self.job_script_filename}')
        import carbonate_job_script_gen as script_gen
        script_gen.reset_config()
        script_gen.set_num_nodes(self.__config_dict['nn'])
        script_gen.set_num_proc_per_node(self.__config_dict['ppn'])
        script_gen.set_job_name(self.__config_dict['jn'])
        script_gen.set_walltime(self.__config_dict['wh'], self.__config_dict['wm'])
        script_gen.set_virtual_mem(self.__config_dict['vmem'])
        script_gen.set_keep_output(keep_output)
        script_gen.set_shell_scripts(self.shell_script_filename)
        if not script_gen.validate_config():
            return None
        self._write_linux(self.job_script_filename, script_gen.run_script())

    def issue_job(self):
        self.generate_input_json()
        self.generate_shell_script()
        self.generate_job_script()

        import subprocess
        import multiprocessing
        p = multiprocessing.Process(target=subprocess.run, args=(["sbatch", self.job_script_filename],))
        p.start()


class SimRunAsyncCarbonate(CoV2VTMSimRunAsync):
    def __init__(self, carbonate_config: dict, root_output_folder, output_frequency=0,
                 screenshot_output_frequency=0, num_runs=1, sim_input=None):
        super().__init__(root_output_folder=root_output_folder, output_frequency=output_frequency,
                         screenshot_output_frequency=screenshot_output_frequency, num_runs=num_runs,
                         sim_input=sim_input)

        self.carbonate_config = carbonate_config

        self.model_config = list()
        for r in range(num_runs):
            model_config = dict()
            model_config['simulation_fname'] = simulation_fname
            model_config['output_frequency'] = output_frequency
            model_config['screenshot_output_frequency'] = screenshot_output_frequency
            model_config['output_dir'] = self.get_run_output_dir(r)
            model_config['sim_input'] = self._sim_input[r]
            self.model_config.append(model_config.copy())

    def run_next_async(self, no_purge=False):
        if not self.has_more_runs:
            return False

        run_idx = self.get_status().index(0, )
        self.run_status[run_idx] = 2

        model_config = self.model_config[run_idx].copy()
        carbonate_config = self.carbonate_config.copy()
        carbonate_config['jn'] += f'_{run_idx}'
        dispatcher = CallableCC3DCarbonateDispatcher(model_config=model_config,
                                                     carb_config=carbonate_config)
        dispatcher.issue_job()
        self.write_sim_inputs(run_idx)
        return True


class CallableCC3DCarbonateScheduler(CallableCoV2VTMScheduler):
    def __init__(self, carbonate_config, root_output_folder, output_frequency=0,
                 screenshot_output_frequency=0, num_runs=1, sim_input=None):
        super().__init__(root_output_folder=root_output_folder, output_frequency=output_frequency,
                         screenshot_output_frequency=screenshot_output_frequency, num_workers=int(1E6),
                         num_runs=num_runs, sim_input=sim_input, dump_dir=None, async_delay=0.001)

        assert isinstance(carbonate_config, dict) or isinstance(carbonate_config, list)
        if isinstance(carbonate_config, dict):
            self.carbonate_config = [carbonate_config.copy() for _ in range(self.num_sets)]
        else:
            self.carbonate_config = carbonate_config
        for x in self.carbonate_config:
            assert isinstance(x, dict)

    def run_instance(self, _set_idx):
        sim_input = self._sim_input[_set_idx].copy()
        BatchRunLib.append_auto_inputs(sim_input)
        return SimRunAsyncCarbonate(carbonate_config=self.carbonate_config[_set_idx],
                                    root_output_folder=self.output_set_directory(_set_idx),
                                    output_frequency=self.output_frequency,
                                    screenshot_output_frequency=self.screenshot_output_frequency,
                                    num_runs=self.num_runs[_set_idx],
                                    sim_input=sim_input)
