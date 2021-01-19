"""\

Usage:
    hprfe2 [-v] [-r PATH] sample [-t PATH] [-a NUM|-s FILE]
    hprfe2 [-v] [-r PATH] sample learn

Arguments:
    -t PATH --template=PATH   Path to a directory with template sampling case files
    -s FILE --strain=FILE     Path to a strain set file
    -a NUM --auto-strain=NUM  Generates strain set of NUM vectors
                              (in the positive quadrant)
    -v                        Verbose output
    -r PATH --root=PATH       Specify the root path of the project, where the
                              configuration file must be located [default: .]

Commands:
    learn                     Run optimization steps

Creates file structure for the sampling. If a path is pass with the -t option, 
it copies template from path, else, it assumes files are already present in the
sampling directory. It requires a file with strain vectors (Voigt notation) for
the generation of the cases. If an integer value between 1 and 63 is passed with 
the -a option, it generates a strain file.

Passing the 'learn' command performs a detection of the elastic range, and adjust 
the timestep accodingly.
"""

import logging
import os
import json
from pathlib import Path
from common import Common


logger = logging.getLogger(__name__)


STRAIN_FN = ["strain_set.dat"]
TEMPL_FN = [
    "MainKratos.py",
    "model.mdpa",
    "materials.json",
    "ProjectParameters.json",
]


def create_case_dir(case, strain, validation=False):
    # create dest dir
    case.mkdir(exist_ok=True)

    # customize properties
    m_prop = case.parent / "ProjectParameters.json"  # template properties file
    p = json.loads(m_prop.read_text())
    p["processes"]["loads_process_list"][0]["Parameters"]["initial_strain"] = strain
    c_prop = case / "ProjectParameters.json"  # destination case properties file
    c_prop.write_text(json.dumps(p, indent=4))
    # customize no-output properties (for speedup calc) in validation cases
    if validation:
        p["processes"]["my_processes"] = []
        p["output_processes"] = {}
        c_prop = case / "ProjectParameters_quiet.json"
        c_prop.write_text(json.dumps(p, indent=4))

    # copy MainKratos.py
    src = case.parent / "MainKratos.py"
    dest = case / "MainKratos.py"
    dest.write_text(src.read_text())

    # copy model.mdpa
    src = case.parent / "model.mdpa"
    dest = case / "model.mdpa"
    dest.unlink(missing_ok=True)  # Remove it before hard-linking it
    src.link_to(dest)  # Create hard link to save space (instead of copy)

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
printf -v ID "%02d" $SLURM_ARRAY_TASK_ID
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


def deploy(common, args):
    # CHECKS:

    # if no template path, files must be present
    if args["--template"] is None:
        for f in TEMPL_FN:
            pf = common.training_path / f
            if not pf.exists():
                logger.error(
                    "No local file '{}'. Missing template path? Aborting.".format(pf)
                )
                exit()
    # if template path, source file set must exist
    else:
        path = Path(args["--template"])
        for f in TEMPL_FN:
            pf = path / f
            if not pf.exists():
                logger.error("Missing remote template file '{}'. Aborting.".format(pf))
                exit()
    # if --strain file passed, file must exist
    n = 0
    if args["--strain"] is not None:
        pf = Path(args["--strain"])
        if not pf.exists():
            logger.error("Missing remote strain file '{}'. Aborting.".format(pf))
            exit()
        n = len(pf.read_text().splitlines())

    # if autostrain option, number of strains is not <1 or >63
    if args["--auto-strain"] is not None:
        s = int(args["--auto-strain"])
        if s < 1 or s > 63:
            logger.error(
                "Selected number of strains ({}) ".format(s)
                + "must be between 1 and 63. Aborting."
            )
            exit()
        if n < s:
            n = s

    # if no remote strain file or autostrain option, strain_set file must exist
    if args["--strain"] is None and args["--auto-strain"] is None:
        pf = common.training_path / STRAIN_FN[0]
        if not pf.exists():
            logger.error(
                "No strain file '{}'present. ".format(pf)
                + "Use -s or -a options to generate one. Aborting."
            )
            exit()
        n = len(pf.read_text().splitlines())

    # number of strains is greater that validation cases
    v = max(common.config["validation_dataset"])
    if n < v:
        logger.error(
            "Number of cases ({}) is smaller than ".format(n)
            + "validation case {}. ".format(v)
            + "Aborting."
        )
        exit()

    # TASKS

    # Create directory if not present already
    if not common.training_path.exists():
        common.training_path.mkdir()
        logger.info("Created directory {}".format(common.training_path))

    # If we have a source path, copy template files
    if args["--template"] is not None:
        path = Path(args["--template"])
        for f in TEMPL_FN:
            src = path / f
            dest = common.training_path / f
            dest.write_text(src.read_text())
        logger.info("Template files copied to sampling directory")
    else:
        logger.info("Using existing template files in directory")

    # If --strain file passed, copy it
    if args["--strain"] is not None:
        src = Path(args["--strain"])
        dest = common.training_path / STRAIN_FN[0]
        dest.write_text(src.read_text())
        logger.info("Strain file copied to sampling directory")

    # If --auto-strain option, create strain set file
    if args["--auto-strain"] is not None:
        n = int(args["--auto-strain"])
        text = ""
        for i in range(1, n + 1):
            b = bin(i)[2:].zfill(6)[::-1]  # binary i, padded with zeroes, and reversed
            line = "{}  {}  {}  {}  {}  {}\n".format(*b)
            text += line
        dest = common.training_path / STRAIN_FN[0]
        dest.write_text(text)
        logger.info("Strain file generated ({} vectors)".format(n))

    # Deploy file structure for sampling
    src = common.training_path / STRAIN_FN[0]
    strain_set = src.read_text().splitlines()
    for i, line in enumerate(strain_set):
        validation = False
        case_path = common.training_path / common.case_name(i)
        strain_vector = [float(x) for x in line.split()]
        if i in common.config["validation_dataset"]:
            validation = True
            logger.info("Case {} set as validation case".format(i))
        create_case_dir(case_path, strain_vector, validation)
        create_run_script(case_path)
        logger.debug("{} {}".format(case_path.name, strain_vector))
    logger.info("Created {} sampling cases".format(i + 1))
    return


def launch_scripts(common):
    # Write launcher scripts
    src = common.training_path / STRAIN_FN[0]
    strain_set = src.read_text().splitlines()
    for i, line in enumerate(strain_set):
        validation = False
        case_path = common.training_path / common.case_name(i)
        strain_vector = [float(x) for x in line.split()]
        if i in common.config["validation_dataset"]:
            validation = True
        create_run_script(case_path)
        logger.debug("{} {}".format(case_path.name, strain_vector))
    create_launchers(common.training_path)
    logger.info("Writen launch scripts".format(i + 1))
    return


def detect_elastic_range(case):
    import KratosMultiphysics
    import KratosMultiphysics.StructuralMechanicsApplication
    import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis as periodic_bc_analysis
    # import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis_xz as periodic_bc_analysis

    param = json.loads((case / "ProjectParameters.json").resolve().read_text())

    #  Replace processes with is_elastic
    param["processes"]["my_processes"] = {}
    block = {
        "Parameters": {
            "model_part_name": "Microstructure.RVE",
            "filename": str(case / "is_elastic.tmp"),
        },
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "IsElastic",
        "python_module": "kratos_process_is_elastic",
    }
    param["processes"]["my_processes"] = [block]
    #  Adapt model and material files
    param["solver_settings"]["model_import_settings"]["input_filename"] = str(
        (case / "model").resolve()
    )
    param["solver_settings"]["material_import_settings"]["materials_filename"] = str(
        (case / "materials.json").resolve()
    )

    # Get params and perform binary search
    t0 = param["problem_data"]["start_time"]
    t1 = param["problem_data"]["end_time"]

    text = ""
    min_dt = (t1 - t0) / 100  # To find it in ~5 iterations
    max_iter = 10
    i = 0
    while (t1 - t0) > min_dt and i < max_iter:
        t = t0 + (t1 - t0) / 2

        # Prepare case
        param["problem_data"]["end_time"] = t
        param["solver_settings"]["time_stepping"] = {}
        param["solver_settings"]["time_stepping"]["time_step"] = t
        parameters = KratosMultiphysics.Parameters(json.dumps(param))
        model = KratosMultiphysics.Model()
        simulation = periodic_bc_analysis.PBCAnalysis(model, parameters)
        simulation.Run()

        elastic = False
        if int((case / "is_elastic.tmp").resolve().read_text()) == 1:
            elastic = True

        if elastic:
            t0 = t
        else:
            t1 = t

        line = "{:.4f} - {:.4f} - Elastic: {}\n".format(t0, t1, elastic)
        i += 1
        text += line
    text += "{:.2f}\n".format(t)
    (case / "elastic.dat").write_text(text)
    (case / "is_elastic.tmp").unlink()


def optimize_timesteps(case, te):

    param = json.loads((case / "ProjectParameters.json").resolve().read_text())

    t0 = param["problem_data"]["start_time"]
    tf = param["problem_data"]["end_time"]
    dt = (tf - t0) / 40  # HARCODED, optimization parameter

    #
    #       fix, 5 ts        fine, dt/2                 coarse, dt*2
    #   t0  .  .  .  . te  . ..t1..... . .   .   t2   .   .   .   .   tf
    #

    ts_table = [
        [t0, (te - t0) / 4],
        [te - (te - t0) * 0.2, (te - t0) / 4],
        [te + (tf - te) * 0.0, dt / 2],
        [te + (tf - te) * 0.4, dt * 2],
    ]

    param["solver_settings"]["time_stepping"] = {}
    param["solver_settings"]["time_stepping"]["automatic_time_step"] = False
    param["solver_settings"]["time_stepping"]["time_step_table"] = ts_table
    (case / "ProjectParameters.json").resolve().write_text(json.dumps(param, indent=4))


def learn(common, args):
    # Deploy file structure for sampling
    src = common.training_path / STRAIN_FN[0]
    strain_set = src.read_text().splitlines()

    # Detect range
    for i, line in enumerate(strain_set):
        case_path = common.training_path / common.case_name(i)
        if (case_path / "elastic.dat").exists():
            logger.info(
                "File {} exists. Skipping case.".format(str(case_path / "elastic.dat"))
            )
            continue
        detect_elastic_range(case_path)

    # Validate results
    for i, line in enumerate(strain_set):
        case_path = common.training_path / common.case_name(i)
        tc = float((case_path / "elastic.dat").read_text().splitlines()[-1])

        if tc < 0.01 or tc > 0.99:
            logger.warning(
                "Elastic region for case {} is outside 0.01 - 0.99 range.".format(
                    str(common.case_name(i))
                )
            )
        print(tc)

    # Adapt params
    for i, line in enumerate(strain_set):
        case_path = common.training_path / common.case_name(i)
        tc = float((case_path / "elastic.dat").read_text().splitlines()[-1])
        optimize_timesteps(case_path, tc)
    logger.info("Optimized timesteps in ProjectParameters.json")

    return


def run(common, args):
    if not args["learn"]:
        deploy(common, args)
        launch_scripts(common)
    else:
        learn(common, args)
