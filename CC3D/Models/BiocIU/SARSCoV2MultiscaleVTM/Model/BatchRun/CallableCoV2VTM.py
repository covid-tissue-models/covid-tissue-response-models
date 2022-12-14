import json
import multiprocessing
import os
import shutil
import sys
import time
sys.path.append(os.environ['PYTHONPATH'])  # Apparently necessary for Linux

from cc3d.CompuCellSetup.CC3DCaller import CC3DCaller, CC3DCallerWorker

sys.path.append(os.path.dirname(os.path.dirname(__file__)))
sys.path.append(os.path.dirname(__file__))

from nCoVToolkit import nCoVUtils
from BatchRunLib import cc3d_input_key

simulation_fname = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'ViralInfectionVTM.cc3d')
simulation_fname = r"/N/u/jferrari/Carbonate/covid_models/CC3D/Models/BiocIU/SARSCoV2MultiscaleVTM/Model/Models/Motion/UniCell_Explicit_Force_2D_with_beta.cc3d"
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

        self._sim_input = sim_input

        # Project convention:   all callable simulation inputs are passed in a dictionary
        #                       key is name of simulation input
        #                       value is value of simulation input
        if isinstance(sim_input, list):
            assert len(sim_input) == num_runs, "Number of runs does not match number of simulation inputs"
        elif isinstance(sim_input, dict):
            print("CoV2VTMSimRun is applying uniform simulation inputs to {} runs".format(self.num_runs))
            self._sim_input = [sim_input] * self.num_runs

        self.sim_output = [None] * self.num_runs

    def set_run_inputs(self, run_idx, sim_inputs):
        assert isinstance(sim_inputs, dict)
        self._sim_input[run_idx] = sim_inputs

    def get_run_output_dir(self, run_idx):
        return os.path.join(self.output_dir_root, f'run_{run_idx}')

    def get_trial_dirs(self):
        return [self.get_run_output_dir(x) for x in range(self.num_runs)]

    def write_sim_inputs(self, run_idx):
        if self._sim_input is None:
            return
        sim_inputs = self._sim_input[run_idx]
        nCoVUtils.export_parameters(sim_inputs, os.path.join(self.get_run_output_dir(run_idx), 'CallableSimInputs.csv'))

    def generate_callable(self, run_idx=0):
        if self._sim_input is not None:
            _sim_input = {cc3d_input_key: self._sim_input[run_idx]}
        else:
            _sim_input = None

        cc3d_caller = CC3DCaller(cc3d_sim_fname=simulation_fname,
                                 output_frequency=self.output_frequency,
                                 screenshot_output_frequency=self.screenshot_output_frequency,
                                 output_dir=self.get_run_output_dir(run_idx),
                                 result_identifier_tag=run_idx,
                                 sim_input=_sim_input)
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
        workers = [CC3DCallerWorker(tasks, results) for _ in range(cov2_vtm_sim_run.num_workers)]
        [w.start() for w in workers]

        # Enqueue jobs
        [tasks.put(cov2_vtm_sim_run.generate_callable(run_idx)) for run_idx in run_list]

        # Add a stop task for each of worker
        [tasks.put(None) for _ in workers]

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


class CoV2VTMSimRunAsync(CoV2VTMSimRun):
    def __init__(self, root_output_folder=generic_root_output_folder, output_frequency=0, screenshot_output_frequency=0,
                 num_runs=1, sim_input=None):
        super().__init__(root_output_folder=root_output_folder, output_frequency=output_frequency,
                         screenshot_output_frequency=screenshot_output_frequency, num_workers=1, num_runs=num_runs,
                         sim_input=sim_input)

        # 0: not yet run; 1: running; 2: done
        self.run_status = multiprocessing.Array('i', [0] * num_runs)

    @property
    def has_more_runs(self) -> bool:
        return any([x == 0 for x in self.run_status])

    @property
    def num_running(self):
        return len([x for x in self.run_status if x == 1])

    @property
    def is_done(self) -> bool:
        return not any([x == 0 or x == 1 for x in self.run_status])

    def set_status(self, _run_idx, _status):
        self.run_status[_run_idx] = _status

    def get_status(self):
        return list(self.run_status)

    def run_next_async(self):
        if not self.has_more_runs:
            return

        run_idx = self.get_status().index(0, )
        self.run_status[run_idx] = 1
        p = multiprocessing.Process(target=sim_run_single, args=(self, run_idx, self.run_status))
        p.start()
        return run_idx


def sim_run_single(cov2_vtm_sim_run: CoV2VTMSimRunAsync, run_idx, status: multiprocessing.Array = None):
    while True:
        # Start worker
        tasks = multiprocessing.JoinableQueue()
        results = multiprocessing.Queue()
        worker = CC3DCallerWorker(tasks, results)
        worker.start()

        # Enqueue jobs
        tasks.put(cov2_vtm_sim_run.generate_callable(run_idx))

        # Add a stop task for each of worker
        tasks.put(None)

        # Monitor worker state
        monitor_rate = 1
        while worker.is_alive():
            time.sleep(monitor_rate)

        # Fetch available results
        while True:
            result = results.get()
            run_idx = result['tag']
            sim_output = result['result']

            print(f'Got CoV2VTMSimRun async result {run_idx}')

            cov2_vtm_sim_run.sim_output[run_idx] = sim_output
            cov2_vtm_sim_run.write_sim_inputs(run_idx)

            if results.empty():
                break

        print(f'CC3DCallerWorker {worker.name} finished with exit code {worker.exitcode}')

        if worker.exitcode == 0:
            if status is not None:
                status[run_idx] = 2
            return


class CallableCoV2VTMScheduler:
    def __init__(self, root_output_folder=generic_root_output_folder, output_frequency=0, screenshot_output_frequency=0,
                 num_workers=1, num_runs=1, sim_input=None, dump_dir=None, async_delay=0):
        self.output_dir_root = root_output_folder
        self.status_file = os.path.join(self.output_dir_root, "batch_status.json")
        if os.path.exists(self.status_file):
            print('Importing from status file: ', self.status_file)
            status_dict = self.load_status(self.status_file)
            output_frequency = status_dict['output_frequency']
            screenshot_output_frequency = status_dict['screenshot_output_frequency']
            num_runs = status_dict['num_runs']
            sim_input = status_dict['sim_input']
            run_status = status_dict['run_status']
            dump_dir = status_dict['dump_dir']
            set_status = status_dict['set_status']
            async_delay = status_dict['async_delay']
        else:
            run_status = None
            set_status = None

        assert output_frequency >= 0
        assert screenshot_output_frequency >= 0
        if isinstance(num_runs, list):
            assert len(num_runs) > 0
        else:
            assert num_runs > 0
        assert num_workers > 0

        self.output_frequency = output_frequency
        self.screenshot_output_frequency = screenshot_output_frequency
        self.num_workers = num_workers
        self.dump_dir = dump_dir
        self.async_delay = async_delay

        # Do version check; simulation inputs via CallableCC3D is an experimental feature as of CompuCell3D v 4.1.0
        def check_callable_cc3d_compat():
            from cc3d.CompuCellSetup import persistent_globals as pg
            assert 'return_object' in dir(pg), "Support for simulation inputs via CallableCC3D not found!"

        if sim_input is not None:
            check_callable_cc3d_compat()

        if sim_input is None:
            self._sim_input = [sim_input]
        elif isinstance(sim_input, list):
            assert not any(not isinstance(x, dict) for x in sim_input)
            self._sim_input = sim_input
        elif isinstance(sim_input, dict):
            self._sim_input = [sim_input]
        else:
            raise ValueError('Incompatible sim inputs')

        self.num_sets = len(self._sim_input)

        if isinstance(num_runs, list):
            assert not any(not isinstance(x, int) for x in num_runs) and len(num_runs) == self.num_sets
            self.num_runs = num_runs
        elif isinstance(num_runs, int):
            self.num_runs = [num_runs] * self.num_sets
        else:
            raise ValueError('Incompatible number of runs')

        if run_status is None:
            run_status = list()
            [run_status.append([0] * x) for x in self.num_runs]
        self.run_status = run_status
        if set_status is None:
            set_status = [0] * self.num_sets
        # 0: working; 1: complete; 2: dumped
        self.set_status = set_status
        self.dump_status()

    def run_instance(self, _set_idx):
        return CoV2VTMSimRunAsync(root_output_folder=self.output_set_directory(_set_idx),
                                  output_frequency=self.output_frequency,
                                  screenshot_output_frequency=self.screenshot_output_frequency,
                                  num_runs=self.num_runs[_set_idx],
                                  sim_input=self._sim_input[_set_idx])

    def dump_status(self):
        status_dict = dict()
        status_dict['output_frequency'] = self.output_frequency
        status_dict['screenshot_output_frequency'] = self.screenshot_output_frequency
        status_dict['num_runs'] = self.num_runs
        status_dict['sim_input'] = self._sim_input
        status_dict['run_status'] = self.run_status
        status_dict['dump_dir'] = self.dump_dir
        status_dict['set_status'] = self.set_status
        status_dict['async_delay'] = self.async_delay
        with open(self.status_file, 'w') as sf:
            json.dump(status_dict, sf)

    @staticmethod
    def load_status(_status_file):
        status_dict = CallableCoV2VTMScheduler.default_status()
        with open(_status_file, 'r') as sf:
            status_dict.update(json.load(sf))
        return status_dict

    @staticmethod
    def default_status():
        status_dict = dict()
        status_dict['output_frequency'] = 0
        status_dict['screenshot_output_frequency'] = 0
        status_dict['num_runs'] = 1
        status_dict['sim_input'] = None
        status_dict['run_status'] = [0]
        status_dict['dump_dir'] = None
        status_dict['set_status'] = [0]
        status_dict['async_delay'] = 0
        return status_dict

    def output_set_directory(self, _set_idx):
        return os.path.join(self.output_dir_root, f"set_{_set_idx}")

    def output_run_directory(self, _set_idx, _run_idx):
        return os.path.join(self.output_set_directory(_set_idx), f"run_{_run_idx}")

    def dump_set_directory(self, _set_idx):
        if self.is_dumping:
            return os.path.join(self.dump_dir, f"set_{_set_idx}")
        else:
            return None

    def dump_run_directory(self, _set_idx, _run_idx):
        if self.is_dumping:
            return os.path.join(self.dump_set_directory(_set_idx), f"run_{_run_idx}")
        else:
            return None

    def prep(self):
        # Create root output directory if necessary
        if not os.path.isdir(self.output_dir_root):
            print('Creating output directory: ', self.output_dir_root)
            os.makedirs(self.output_dir_root)

        # Create set subdirectories if necessary
        for s in range(self.num_sets):
            if self.set_status[s] == 0 and not os.path.isdir(self.output_set_directory(s)):
                os.mkdir(self.output_set_directory(s))

        # Purge unfinished runs
        for s in range(self.num_sets):
            if self.set_status[s] != 0:
                continue
            for r in range(self.num_runs[s]):
                r_status = self.run_status[s][r]
                r_dir = self.output_run_directory(s, r)
                if os.path.isdir(r_dir) and r_status != 2:
                    shutil.rmtree(r_dir)
                    self.run_status[s][r] = 0

        # Create dump directory if necessary
        if self.is_dumping and not os.path.isdir(self.dump_dir):
            print('Creating dump directory: ', self.dump_dir)
            os.makedirs(self.dump_dir)

    @property
    def fin_key(self) -> int:
        if self.dump_dir is not None:
            return 2
        else:
            return 1

    @property
    def is_dumping(self) -> bool:
        return self.fin_key == 2

    def run(self):
        self.prep()
        set_runners = [self.run_instance(s) for s in range(self.num_sets)]
        for s in range(self.num_sets):
            for r in range(len(self.run_status[s])):
                set_runners[s].set_status(r, self.run_status[s][r])

        # Run block
        import time
        while True:
            if not any([s != self.fin_key for s in self.set_status]):
                if self.is_dumping:
                    self.check_dumps()
                self.dump_status()
                return

            # Check status of every set runner
            for s in range(self.num_sets):
                if self.set_status[s] != 0:
                    if self.is_dumping and self.set_status[s] == 1 and not os.path.isdir(self.output_set_directory(s)):
                        self.set_status[s] = 2
                    continue

                set_runner = set_runners[s]
                if set_runner.is_done:
                    self.set_status[s] = 1
                    if self.is_dumping:
                        move_dir_async(self.output_set_directory(s), self.dump_set_directory(s))

            num_running = sum([sr.num_running for sr in set_runners])
            num_to_add = self.num_workers - num_running
            num_added = 0
            for s in range(self.num_sets):
                while set_runners[s].has_more_runs and num_added < num_to_add:
                    if set_runners[s].run_next_async():
                        num_added += 1
                        if self.async_delay > 0:
                            time.sleep(self.async_delay)
                self.run_status[s] = set_runners[s].get_status()

            # Update status file and then wait a little bit
            self.dump_status()
            time.sleep(10)

    def check_dumps(self):
        # Hold until dumping is complete
        if not self.is_dumping:
            return

        import time
        while any([self.set_status[s] != 2 for s in range(self.num_sets)]):
            for s in range(self.num_sets):
                if self.set_status[s] == 1 and not os.path.isdir(self.output_set_directory(s)):
                    self.set_status[s] = 2
                    self.dump_status()
            time.sleep(10)


class _MoveDirProcess(multiprocessing.Process):
    def __init__(self, _src_dir, _tgt_dir):
        """
        Process to asynchronously move one directory to another
        :param _src_dir: directory to move
        :param _tgt_dir: location of resulting move
        """
        super().__init__()
        self._src_dir = _src_dir
        self._tgt_dir = _tgt_dir

    def run(self):
        shutil.move(self._src_dir, self._tgt_dir)


def move_dir_async(_src_dir, _tgt_dir) -> None:
    p = _MoveDirProcess(_src_dir, _tgt_dir)
    p.start()
