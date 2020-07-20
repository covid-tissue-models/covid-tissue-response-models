import multiprocessing
import os
import time

from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller, CC3DCallerWorker

from nCoVToolkit import nCoVUtils
from BatchPostCoV2VTM import CallableCC3DRenderer, CoV2VTMSimRunPost

simulation_fname = os.path.join(os.path.dirname(__file__), 'ViralInfectionVTM.cc3d')
generic_root_output_folder = os.path.abspath(os.path.join(os.path.splitdrive(os.getcwd())[0], '/CallableCoV2VTM'))


class CoV2VTMSimRun:
    def __init__(self, root_output_folder=generic_root_output_folder, output_frequency=0, screenshot_output_frequency=0,
                 num_workers=1, num_runs=1, sim_input=None):

        assert output_frequency >= 0
        assert screenshot_output_frequency >= 0
        assert num_runs > 0
        assert num_workers > 0

        self.output_frequency = output_frequency
        self.screenshot_output_frequency = screenshot_output_frequency
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
        return os.path.join(self.output_dir_root, f'run_{run_idx}')

    def get_trial_dirs(self):
        return [self.get_run_output_dir(x) for x in range(self.num_runs)]

    def write_sim_inputs(self, run_idx):
        if self.__sim_input is None:
            return
        sim_inputs = self.__sim_input[run_idx]
        nCoVUtils.export_parameters(sim_inputs, os.path.join(self.get_run_output_dir(run_idx), 'CallableSimInputs.csv'))

    def generate_callable(self, run_idx=0):
        if self.__sim_input is not None:
            sim_input = self.__sim_input[run_idx]
        else:
            sim_input = None

        cc3d_caller = CC3DCaller(cc3d_sim_fname=simulation_fname,
                                 output_frequency=self.output_frequency,
                                 screenshot_output_frequency=self.screenshot_output_frequency,
                                 output_dir=self.get_run_output_dir(run_idx),
                                 result_identifier_tag=run_idx,
                                 sim_input=sim_input)
        return cc3d_caller


def run_cov2_vtm_sims(cov2_vtm_sim_run: CoV2VTMSimRun) -> CoV2VTMSimRun:
    # Make complete list of jobs
    run_list = [run_idx for run_idx in range(cov2_vtm_sim_run.num_runs)]

    while True:
        num_jobs_curr = len(run_list)
        print('Doing CoV2VTMSimRun batch iteration with {} remaining jobs.'.format(num_jobs_curr))

        # Start workers
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        workers = [CC3DCallerWorker(tasks, results) for i in range(cov2_vtm_sim_run.num_workers)]
        [w.start() for w in workers]

        # Enqueue jobs
        [tasks.put(cov2_vtm_sim_run.generate_callable(run_idx)) for run_idx in run_list]

        # Add a stop task for each of worker
        [tasks.put(None) for w in workers]

        # Monitor worker state
        monitor_rate = 1
        while [w for w in workers if w.is_alive()]:
            time.sleep(monitor_rate)

        # Fetch available results
        while True:
            result = results.get()
            run_idx = result['tag']
            sim_output = result['result']

            print('Got CoV2VTMSimRun batch result {}'.format(run_idx))

            cov2_vtm_sim_run.sim_output[run_idx] = sim_output
            cov2_vtm_sim_run.write_sim_inputs(run_idx)
            run_list.remove(run_idx)

            if results.empty():
                break

        num_jobs_next = len(run_list)
        print('CoV2VTMSimRun batch results processed with {} remaining jobs.'.format(len(run_list)))
        [print('CC3DCallerWorker {} finished with exit code {}'.format(w.name, w.exitcode)) for w in workers]

        if not run_list:
            print('CoV2VTMSimRun batch complete!')
            break
        elif num_jobs_next == num_jobs_curr:
            print('CoV2VTMSimRun batch run failed! Terminating early.')
            break

    return cov2_vtm_sim_run


# Example of usage / convenience sequence to do intended overall workflow
# if __name__ == '__main__':
#     # Setup batch run
#     _root_output_folder = generic_root_output_folder
#     cov2_vtm_sim_run = CoV2VTMSimRun(num_runs=10,
#                                      num_workers=5,
#                                      output_frequency=10,
#                                      screenshot_output_frequency=10,
#                                      root_output_folder=_root_output_folder)
#
#     # Execute batch simulations
#     cov2_vtm_sim_run = run_cov2_vtm_sims(cov2_vtm_sim_run)
#
#     # Export model parameters
#     from Simulation import ViralInfectionVTMModelInputs
#     export_file_abs = os.path.join(os.path.abspath(cov2_vtm_sim_run.output_dir_root),
#                                    "ViralInfectionVTMModelParams.csv")
#     nCoVUtils.export_parameters(ViralInfectionVTMModelInputs, export_file_abs)
#
#     # Post-process metrics
#     cov2_vtm_sim_run_post = CoV2VTMSimRunPost(cov2_vtm_sim_run)
#     cov2_vtm_sim_run_post.export_transient_plot_trials()
#     cov2_vtm_sim_run_post.export_transient_plot_stat()
#     cov2_vtm_sim_run_post.export_2var_plot_trials('MedViral', 'MedCyt')
#     cov2_vtm_sim_run_post.export_2var_plot_stat('MedViral', 'MedCyt', plot_stdev=False)
#
#     # Render field data
#     callable_cc3d_renderer = CallableCC3DRenderer(cov2_vtm_sim_run)
#     screenshot_specs = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'screenshots.json')
#     callable_cc3d_renderer.load_screenshot_specs(screenshot_specs)
#     callable_cc3d_renderer.render_results()
