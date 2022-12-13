import os

DEBUG = "debug"
INTERACTIVE = "interactive"


# Per general compute nodes
_ciu_lims = {'nn': 72,
             'ppn': 24,
             'vmem': 503,
             'wh': None}

# Debug queue
_debug_lims = {'nn': 2,
               'ppn': 24,
               'vmem': 503,
               'wh': 1}

# Interactive queue
_interact_lims = {'nn': 1,
                  'ppn': 8,
                  'vmem': 503,
                  'wh': 8}


_default_config = {'ec': None,
                   'ee': None,
                   'jn': "MYJOB",
                   'ko': True,
                   'nn': 1,
                   'ppn': 1,
                   'shell_scripts': None,
                   'vmem': None,
                   'wh': None,
                   'wm': None,
                   'q': None}
_config = {k: v for k, v in _default_config.items()}
_hw_lims = {k: v for k, v in _ciu_lims.items()}


def reset_config_entry(_k: str):
    _config[_k] = _default_config[_k]


def reset_config():
    for _k in _default_config.keys():
        _config[_k] = _default_config[_k]


def set_email_contact(_ec: str):
    _config['ec'] = _ec


def reset_email_contact():
    reset_config_entry('ec')


def set_email_events(_ee: str):
    check_events(_ee)
    _config['ee'] = _ee


def reset_email_events():
    reset_config_entry('ee')


def set_job_name(_jn: str):
    _config['jn'] = _jn


def reset_job_name():
    reset_config_entry('jn')


def set_keep_output(_ko: bool):
    _config['ko'] = _ko


def reset_keep_output():
    reset_config_entry('ko')


def set_num_nodes(_nn):
    _config['nn'] = _nn


def reset_num_nodes():
    reset_config_entry('nn')


def set_num_proc_per_node(_ppn):
    _config['ppn'] = _ppn


def reset_num_proc_per_node():
    reset_config_entry('ppn')


def set_shell_scripts(_ss):
    if isinstance(_ss, str):
        add_shell_script(_ss)
    elif isinstance(_ss, list):
        [add_shell_script(_s) for _s in _ss]
    else:
        raise ValueError(f"Invalid shell script input. Can be a string or list of strings")


def add_shell_script(_ss: str):
    if not isinstance(_ss, str):
        raise ValueError(f"Invalid shell script input: {_ss}")
    if _config['shell_scripts'] is None:
        _config['shell_scripts'] = list()
    _config['shell_scripts'].append(_ss)


def reset_shell_scripts():
    reset_config_entry('shell_scripts')


def set_virtual_mem(_vmem):
    _config['vmem'] = _vmem


def reset_virtual_mem():
    reset_config_entry('vmem')


def set_walltime(hours: int = None, minutes: int = None):
    if hours is None and minutes is None:
        return
    _config['wh'] = hours
    _config['wm'] = minutes


def reset_walltime():
    reset_config_entry('wh')
    reset_config_entry('wm')


def check_events(_ee: str):
    if not len(_ee.replace("a", "").replace("b", "").replace("e", "")) == 0:
        raise ValueError("Contact events must be specified with a, b and/or e")


def set_queue(debug: bool = False, interactive: bool = False):
    if debug:
        _hw_lims.update(_debug_lims)
        _config['q'] = DEBUG
    elif interactive:
        _hw_lims.update(_interact_lims)
        _config['q'] = INTERACTIVE


def optimize():
    if not any([_config[k] > _interact_lims[k] for k in _interact_lims.keys()]):
        set_queue(interactive=True)
    elif not any([_config[k] > _debug_lims[k] for k in _debug_lims.keys()]):
        set_queue(debug=True)


def keep_job_output(_ko):
    if _ko:
        return "#PBS -k o\n"
    else:
        return ""


def job_name(_jn):
    assert len(_jn) > 0, "Job needs a name"
    return f"#SBATCH -J {_jn}\n"


def email_contact(_ec):
    return f"#PBS -M {_ec}\n"


def job_requirements(_nn, _ppn, vmem=None, hours=None, minutes=None):
    assert _nn > 0
    assert _ppn > 0
    # o = f"#PBS -l nodes={_nn}:ppn={_ppn}"
    o = f"#SBATCH  --nodes={_nn}\n"
    o += f"#SBATCH --ntasks-per-node={_ppn}\n"
    if walltime(hours, minutes) is not None:
        o += f"#SBATCH {walltime(hours, minutes)}\n"
    if virtual_mem(vmem) is not None:
        o += f"#SBATCH {virtual_mem(vmem)}"
    o += "\n"
    return o


def queue(_q: str):
    return "\n"
    if _q == DEBUG:
        return "#PBS -q " + DEBUG + "\n"
    else:
        return "\n"


def walltime(hours: int = None, minutes: int = None):
    if hours is not None or minutes is not None:
        if hours is None:
            print(f"Carbonate jobs are allotted 60 minutes of wall time by default. Ignoring input: {minutes} minutes")

        if minutes is None:
            minutes = "00"
        elif minutes < 10:
            minutes = f"0{minutes}"
        else:
            minutes = str(minutes)
        if hours < 10:
            hours = f"0{hours}"
        else:
            hours = str(hours)
        return f"--time={hours}:{minutes}:00"
    return None


def virtual_mem(_vmem):
    if _vmem is not None:
        return f"--mem={_vmem}G"
    return None


def email_events(email_event_labs: str = None, email_addr: str = None):
    if email_event_labs is not None and email_addr is not None:
        check_events(email_event_labs)
        return f"#PBS -M {email_addr}\n" + f"#PBS -m {email_event_labs}\n"
    else:
        return ""


def join_outputs():
    return "#PBS -j oe\n"


def targets(_shell_scripts):
    if isinstance(_shell_scripts, list):
        if len(_shell_scripts) > 1:
            _str = ""
            for t in _shell_scripts:
                _str += f"sh {t} &\n"
        else:
            _str = "sh " + _shell_scripts[0] + "\n"
        return _str
    elif isinstance(_shell_scripts, str):
        return _shell_scripts + "\n"
    assert False, 'No valid target scripts specified'


def validate_config() -> bool:
    """
    Enforces Carbonate hardware constraints
    :return: True if valid configuration
    """
    if _config['nn'] > _hw_lims['nn']:
        print(f"Number of nodes exceeds Carbonate limits ({_hw_lims['nn']})")
        return False
    if _config['ppn'] > _hw_lims['ppn']:
        print(f"Processors per node exceeds Carbonate limits ({_hw_lims['ppn']})")
        return False
    if _config['vmem'] is not None and _config['vmem'] > _hw_lims['vmem']:
        print(f"Virtual memory exceeds Carbonate limits ({_hw_lims['vmem']})")
        return False
    if _hw_lims['wh'] is not None:
        if _config['wh'] > _hw_lims['wh'] or \
                (_config['wh'] == _hw_lims['wh'] and _config['wm'] is not None and _config['wm'] > 0):
            print(f"Walltime exceeds Carbonate {_debug_lims['q']} queue limits ({_hw_lims['wh']})")
            return False
    return True


def interactive_run_script() -> str:
    cmd = "qsub -I -q interactive -l "
    cmd += f"nodes={_config['nn']}:ppn={_config['ppn']},{walltime(_config['wh'], _config['wm'])}"
    o = f"CMD=\"{cmd}\"\n"
    o += "eval $CMD\n"
    o += targets(_config['shell_scripts'])
    return o

def out_file(filename):
    s = f"# SBATCH -o {filename}_%j.txt\n"
    s += f"# SBATCH -e {filename}_%j.err\n"
    return s


def run_script(_jo: bool = True, job_queue: str = "", opt_queue: bool = False) -> str:
    if opt_queue:
        optimize()
    if job_queue == DEBUG or job_queue == INTERACTIVE:
        set_queue(debug=job_queue == DEBUG, interactive=job_queue == INTERACTIVE)
    if _config['q'] == INTERACTIVE:
        return interactive_run_script()

    o = "#!/bin/bash\n"
    # o += keep_job_output(_config['ko'])
    o += job_requirements(_config['nn'], _config['ppn'], _config['vmem'], _config['wh'], _config['wm'])
    # o += queue(_config['q'])
    # o += email_events(_config['ee'], _config['ec'])
    o += job_name(_config['jn'])
    o += out_file(_config['jn'])
    # if _jo:
    #     o += join_outputs()
    o += targets(_config['shell_scripts'])
    return o


def readme() -> str:
    return """
Per https://kb.iu.edu/d/aolp,  

"
Carbonate has 72 general-purpose compute nodes, each with 256 GB of RAM, and eight large-memory compute nodes, each 
with 512 GB of RAM. Each general-purpose compute node is a Lenovo NeXtScale nx360 M5 server equipped with two 12-core 
Intel Xeon E5-2680 v3 CPUs and four 480 GB solid-state drives.
"

The following is an example of using this script generator

#           Job details

# List or single entry of shell scripts to execute in this job
#   Paths are relative to directory from where this script will be called, or absolute
shell_scripts = [r"/geode2/home/u100/tjsego/Carbonate/batch_run_carbonate.sh"]

#           Job execution details

# Number of nodes
nn = 2
# Number of processors per node
ppn = 10
# Job name
jn = "FREEZEHSOLO"
# Email contact (optional, None)
ec = r"jthutt@tatooine.net"
# Contact events (a = aborted, b = begin, e = end); e.g., abe
ee = "abe"
# Keep job output?
ko = False
# Wall clock time (hours, minutes); necessary if job needs more than 60 minutes
wh = 2
wm = None
# Set to a value in GB if job needs more than 16 GB
vmem = None

import os
import carbonate_job_script_gen as cjsg
cjsg.set_num_nodes(nn)
cjsg.set_num_proc_per_node(ppn)
cjsg.set_job_name(jn)
cjsg.set_email_contact(ec)
cjsg.set_email_events(ee)
cjsg.set_keep_output(ko)
cjsg.set_walltime(wh, wm)
cjsg.set_virtual_mem(vmem)
cjsg.set_shell_scripts(shell_scripts)
cjsg.optimize()
run_script = cjsg.run_script()
fn = os.path.join(os.getcwd(), "job.script")
with open(fn, "w") as f:
    f.write(run_script)
"""


def main():
    print(readme())


if __name__ == "__main__":
    main()
