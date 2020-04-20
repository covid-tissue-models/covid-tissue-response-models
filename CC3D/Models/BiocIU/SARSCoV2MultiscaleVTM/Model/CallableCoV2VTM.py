import csv
import multiprocessing
from os.path import dirname, join
from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller, CC3DCallerWorker

from nCoVToolkit import nCoVUtils

simulation_fname = join(dirname(__file__), 'Coronavirus.cc3d')
generic_root_output_folder = join(dirname(__file__), 'CallableCoV2VTM')


class CoV2VTMSimRun:
    def __init__(self, root_output_folder=generic_root_output_folder, output_frequency=0, screensht_output_frequency=0,
                 num_workers=1, num_runs=1, sim_input=None):

        assert output_frequency >= 0
        assert screensht_output_frequency >= 0
        assert num_runs > 0
        assert num_workers > 0

        self.output_frequency = output_frequency
        self.screensht_output_frequency = screensht_output_frequency
        self.output_dir_root = root_output_folder
        self.num_workers = num_workers
        self.num_runs = num_runs

        # Do version check; simulation inputs via CallableCC3D is an experimental feature as of CompuCell3D v 4.1.0
        def check_callable_cc3d_compat():
            from cc3d.CompuCellSetup import persistent_globals as pg
            assert 'return_object' in dir(pg), "Support for simulation inputs via CallableCC3D not found!"

        if sim_input is not None:
            check_callable_cc3d_compat()

        self.__sim_input = sim_input

        # Project convention:   all callable simulation inputs are passed in a dictionary
        #                       key is name of simulation input
        #                       value is value of simulation input
        if isinstance(sim_input, list):
            assert len(sim_input) == num_runs, "Number of runs does not match number of simulation inputs"
        elif isinstance(sim_input, dict):
            print("CoV2VTMSimRun is applying uniform simulation inputs to {} runs".format(self.num_runs))
            self.__sim_input = [sim_input] * self.num_runs

        self.sim_output = [None] * self.num_runs

    def set_run_inputs(self, run_idx, sim_inputs):
        assert isinstance(sim_inputs, dict)
        self.__sim_input[run_idx] = sim_inputs

    def get_run_output_dir(self, run_idx):
        return join(self.output_dir_root, f'run_{run_idx}')

    def write_sim_inputs(self, run_idx):
        if self.__sim_input is None:
            return
        sim_inputs = self.__sim_input[run_idx]
        nCoVUtils.export_parameters(sim_inputs, join(self.get_run_output_dir(run_idx), 'CallableSimInputs.csv'))

    def generate_callable(self, run_idx=0):
        if self.__sim_input is not None:
            sim_input = self.__sim_input[run_idx]
        else:
            sim_input = None

        cc3d_caller = CC3DCaller(cc3d_sim_fname=simulation_fname,
                                 output_frequency=self.output_frequency,
                                 screenshot_output_frequency=self.screensht_output_frequency,
                                 output_dir=self.get_run_output_dir(run_idx),
                                 result_identifier_tag=run_idx,
                                 sim_input=sim_input)
        return cc3d_caller


def run_cov2_vtm_sims(cov2_vtm_sim_run: CoV2VTMSimRun) -> CoV2VTMSimRun:

    # Start workers
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()
    workers = [CC3DCallerWorker(tasks, results) for i in range(cov2_vtm_sim_run.num_workers)]
    [w.start() for w in workers]

    # Enqueue jobs
    for run_idx in range(cov2_vtm_sim_run.num_runs):
        cc3d_callable = cov2_vtm_sim_run.generate_callable(run_idx)
        tasks.put(cc3d_callable)

    # Add a stop task for each of worker
    for i in range(cov2_vtm_sim_run.num_workers):
        tasks.put(None)

    # Wait for all of the tasks to finish
    tasks.join()

    # Start printing results
    while cov2_vtm_sim_run.num_runs:
        result = results.get()
        run_idx = result['tag']
        sim_output = result['result']
        cov2_vtm_sim_run.sim_output[run_idx] = sim_output
        cov2_vtm_sim_run.write_sim_inputs(run_idx)
        cov2_vtm_sim_run.num_runs -= 1
        print('{} runs remaining'.format(cov2_vtm_sim_run.num_runs))

    return cov2_vtm_sim_run


if __name__ == '__main__':
    cov2_vtm_sim_run = CoV2VTMSimRun(num_runs=10,
                                     num_workers=5,
                                     output_frequency=10,
                                     screensht_output_frequency=10)
    run_cov2_vtm_sims(cov2_vtm_sim_run)
