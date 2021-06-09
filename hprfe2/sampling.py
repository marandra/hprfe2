"""\

Usage:
    hprfe2 [-v] [-r PATH] deploy [-f]

Arguments:
    -v --verbose              Verbose output
    -r PATH --root=PATH       Specify the root path of the project, where the
                              configuration file must be located [default: .]
    -f --fromdb               (Over)write template files from db

Commands:
    deploy                    Create sampling file structure and launch scripts, using
                              existing or provided template files

Creates a file structure for sampling in the 'ROOT_PATH/sampling/' directory.
Template files must be already present or use the -f option to write or overwrite
them with the copy in the database. The template files in the database are updated
every time the files are deployed. Template files require a file with a list of
strain vectors (Voigt notation) for the definition of the training trajectories.
"""

import json
import logging
import os
from pathlib import Path

from common import Common

logger = logging.getLogger(__name__)


TEMPL_FN = {
    "MAIN": "MainKratos.py",
    "MODEL": "model.mdpa",
    "MATERIALS": "materials.json",
    "PARAMETERS": "ProjectParameters.json",
    "STRAINSET": "strain_set.dat",
}


def create_launchers(path):
    """Create auxiliary files for running cases"""

    fname = "launcher_slurm.bash"
    script = """\
#!/bin/bash
#SBATCH --job-name=fiber_array
#SBATCH --ntasks-per-core=1
#SBATCH --ntasks=1
#SBATCH --array=00-40

#SBATCH --partition=HM
##SBATCH --mem-per-cpu=1024
##SBATCH --time=03:00:00

export OMP_NUM_THREADS=1
printf -v ID "%03d" $SLURM_ARRAY_TASK_ID
bash tmp_case_${ID}.bash
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

for SCRIPT in tmp_*.bash
do
  bash $SCRIPT
  sleep 0.05
done
"""
    (path / fname).write_text(script)


class Sampling:
    def __init__(self, common, args):
        self.common = common
        self.args = args

    def check_template(self):
        # Template files must be present
        for f in TEMPL_FN.values():
            pf = self.common.training_path / f
            if not pf.exists():
                logger.error(f"Missing file '{pf}'. Aborting.")
                exit()

        # Number of strains is greater that validation cases
        pf = self.common.training_path / TEMPL_FN["STRAINSET"]
        n = len(pf.read_text().splitlines())
        v = max(self.common.config["validation_dataset"])
        if n < v:
            logger.error(
                f"Number of cases ({n}) is smaller than "
                + f"validation case {v}. "
                + "Aborting."
            )
            exit()

    def save_template(self):
        for k, v in TEMPL_FN.items():
            path = self.common.training_path / v
            data = path.read_text()
            self.common.set_dataset(data, "TEMPLATE", k, replace=True)

    def load_template(self):
        self.common.training_path.mkdir(exist_ok=True)
        for k, v in TEMPL_FN.items():
            data = self.common.get_dataset("TEMPLATE", k)
            path = self.common.training_path / v
            path.write_text(data)

    def generate_cases(self):
        # Generate cases
        src = self.common.training_path / TEMPL_FN["STRAINSET"]
        strain_set = src.read_text().splitlines()
        self.cases = []
        for i, line in enumerate(strain_set):
            strain = [float(x) for x in line.split()]
            case = Case(self.common, i, strain)
            if i in self.common.config["validation_dataset"]:
                case.is_validation = True
                logger.info(f"Case {i} set as validation case")
            self.cases.append(case)

    def deploy_cases(self):
        # Deploy file structure for sampling
        for c in self.cases:
            c.create_directory(self.common)
            c.create_script()
        logger.info(f"Created {len(self.cases)} sampling cases")
        create_launchers(self.common.training_path)
        logger.info("Written launch scripts")


class Case:
    def __init__(self, common, i, strain_vector, is_validation=False):
        self.name = common.case_name(i)
        self.path = common.training_path / self.name
        self.strain = strain_vector
        self.is_validation = is_validation

    def create_directory(self, common):
        # create dest dir
        self.path.mkdir(exist_ok=True)

        # customize properties
        m_prop = self.path.parent / "ProjectParameters.json"  # template properties file
        p = json.loads(m_prop.read_text())
        vect = [f"{v}" for v in self.strain]
        p["processes"]["loads_process_list"][0]["Parameters"]["imposed_strain"] = vect
        # TODO: Fix the path, it need root_path. This is a workaround
        p["processes"]["my_processes"][1]["Parameters"]["material_root_path"] = str(
            common.root_path
        )
        c_prop = (
            self.path / "ProjectParameters.json"
        )  # destination case properties file
        c_prop.write_text(json.dumps(p, indent=4))
        # customize no-output properties (for speedup calc) in validation cases
        if self.is_validation:
            p["processes"]["my_processes"] = []
            p["output_processes"] = {}
            c_prop = self.path / "ProjectParameters_quiet.json"
            c_prop.write_text(json.dumps(p, indent=4))

        # copy MainKratos.py
        src = self.path.parent / "MainKratos.py"
        dest = self.path / "MainKratos.py"
        dest.write_text(src.read_text())

        # link model.mdpa
        src = self.path.parent / "model.mdpa"
        dest = self.path / "model.mdpa"
        dest.unlink(missing_ok=True)  # Remove it before hard-linking it
        src.link_to(dest)  # Create hard link to save space (instead of copy)

        # copy materials.json
        src = self.path.parent / "materials.json"
        dest = self.path / "materials.json"
        dest.write_text(src.read_text())
        # and link rve data if present
        srcs = self.path.parent.glob(common.rve_fname("*", "*", "*"))
        for src in srcs:
            dest = self.path / src.name
            dest.unlink(missing_ok=True)  # Remove it before hard-linking it
            src.link_to(dest)  # Create hard link to save space (instead of copy)

    def create_script(self):
        """
        Writes temporary launch script for each case (to be run externally)
        """
        script_fname = "tmp_" + self.name + ".bash"
        script = ""
        script += "export OMP_NUM_THREADS=1\n"
        script += "export PYTHONPATH={}\n".format(os.environ["PYTHONPATH"])
        script += "export LD_LIBRARY_PATH={}\n".format(os.environ["LD_LIBRARY_PATH"])
        script += "cd {}\n".format(self.name)
        if False:
            script += "learn.py\n"
        script += "python MainKratos.py > outMainKratos\n"
        if self.is_validation:
            script += "/usr/bin/time -v -o time.dat python MainKratos.py "
            script += "ProjectParameters_quiet.json > outMainKratos_quiet\n"
        script += "cd ..\n"
        script += "rm {}\n".format(script_fname)
        (self.path.parent / script_fname).write_text(script)


# def print_histogram(values, bins=20, t0=0, t1=1):
#    text = "Elastic range distribution:\n\n"
#    B = [0] * bins
#    for v in values:
#        b = int(v / (t1 - t0) * bins)
#        B[b] += 1
#    for i, c in enumerate(B):
#        text += "{:0.2f}: {:>3d} |{}\n".format(i * (t1 - t0) / bins, c, "*" * c)
#    return text


def run(common, args):
    sampling = Sampling(common, args)
    if args["--fromdb"]:
        sampling.load_template()
    else:
        sampling.save_template()
    sampling.check_template()
    sampling.generate_cases()
    sampling.deploy_cases()
