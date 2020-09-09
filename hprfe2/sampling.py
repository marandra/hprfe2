"""
docstrings here
"""

import os
import json
from pathlib import Path
from common import Common


def customize_properties(props, strain, quiet=False):
    """
    Updates template properties with strain for case.
    Optionally removes all output (for time measuring).
    Receives properties dict, strain to write, and output flag.
    """

    props["processes"]["loads_process_list"][0]["Parameters"]["initial_strain"] = strain

    if quiet:
        props["processes"]["my_processes"] = []
        props["output_processes"] = {}

    return props


def create_case_dir(case, strain):
    """
    Files and dirs structure:
    case: root_path/training/trajectory_35/
    source files: root_path/training/MainKratos.py
                                     macro_model.mdpa
                                     ProjectParameters.json
    """

    # create dest dir
    case.mkdir(exist_ok=True)

    # customize properties
    m_prop = case.parent / "ProjectParameters.json"  # template properties file
    p = json.loads(m_prop.read_text())
    p = customize_properties(p, strain)
    c_prop = case / "ProjectParameters.json"  # destination case properties file
    c_prop.write_text(json.dumps(p, indent=4))

    # customize no-output properties (for speedup calc)
    p = json.loads(m_prop.read_text())
    p = customize_properties(p, strain, quiet=True)
    c_prop = case / "ProjectParameters_quiet.json"
    c_prop.write_text(json.dumps(p, indent=4))

    # copy MainKratos.py
    src = case.parent / "MainKratos.py"
    dest = case / "MainKratos.py"
    dest.write_text(src.read_text())

    # copy model.mdpa
    src = case.parent / "model.mdpa"
    dest = case / "model.mdpa"
    dest.write_text(src.read_text())

    # copy materials.json
    src = case.parent / "materials.json"
    dest = case / "materials.json"
    dest.write_text(src.read_text())

    return


def create_run_script(case):
    """
    Writes temporary launch script for each case (run externally)
    """

    script_fname = "tmp_" + case.name + ".bash"
    script = """\
export OMP_NUM_THREADS=1
export PYTHONPATH={}
export LD_LIBRARY_PATH={}
cd {}
/usr/bin/time -v -o time.dat python MainKratos.py > outMainKratos
#/usr/bin/time -v -o time_quiet.dat python MainKratos.py ProjectParameters_quiet.json > outMainKratos_quiet
cd ..
rm {}
""".format(
        os.environ["PYTHONPATH"],
        os.environ["LD_LIBRARY_PATH"],
        case.name,
        script_fname,
    )
    (case.parent / script_fname).write_text(script)


def create_launchers(path):
    """Create auxiliary files for running cases"""

    fname = "launcher_slurm.bash"
    script = """\
#!/bin/bash
#SBATCH --job-name=fiber_array
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=1
#SBATCH --array=00-40

## Settings for "fiber2" case:
##
#SBATCH --partition=HM
#SBATCH --mem-per-cpu=1024
#SBATCH --time=03:00:00

## Settings for "fiber3" case:
##
##SBATCH --partition=HM
##SBATCH --mem-per-cpu=????
##SBATCH --time=????

export OMP_NUM_THREADS=1
printf -v ID "%02d\n" $SLURM_ARRAY_TASK_ID
TRAJECTORYPATH=$PWD/case_$ID
cd $TRAJECTORYPATH
# TODO: update the command (output time.dat and time_quiet.dat)
time python3 MainKratos.py
cd ..
mv slurm-$SLURM_ARRAY_JOB_ID\_$SLURM_ARRAY_TASK_ID.out $TRAJECTORYPATH
"""
    (path / fname).write_text(script)

    fname = "launcher_pueue.bash"
    script = """\
#!/bin/bash

for SCRIPT in tmp_*.bash
do
  pueue add -- bash $SCRIPT
  sleep 0.05
done
"""
    (path / fname).write_text(script)

    fname = "launcher.bash"
    script = """\
#!/bin/bash

(for SCRIPT in tmp_*.bash
do
  bash $SCRIPT
  sleep 0.05
done
"""
    (path / fname).write_text(script)


#######################################
# main
#######################################

if __name__ == "__main__":

    import sys

    if len(sys.argv) > 1:
        co = Common(root_path=Path(sys.argv[1]))
    else:
        exit("Missing root_path argument.")

    strain_set = co.parse_training_strain_set()
    for i, line in enumerate(strain_set):
        case_path = co.training_path / co.case_name(i)
        strain_vector = [float(x) for x in line.split()]
        create_case_dir(case_path, strain_vector)
        create_run_script(case_path)
        print(case_path.name, strain_vector)
    create_launchers(co.training_path)
