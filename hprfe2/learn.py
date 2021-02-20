#! /usr/bin/python
import os
import json
from pathlib import Path

import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis as analysis
# import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis_xz as analysis


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
        # Copy RVE materials, if present
        srcs = path.glob(common.rve_fname("*", "*", "*"))
        for src in srcs:
            dest = common.training_path / src.name
            dest.write_bytes(src.read_bytes())
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
        create_case_dir(common, case_path, strain_vector, validation)
        create_run_script(case_path)
        logger.debug("{} {}".format(case_path.name, strain_vector))
    logger.info("Created {} sampling cases".format(i + 1))
    return


def adjust_elastic_range(path, outfile="elastic.dat"):
    """Adjust the time range to fit the elastic region in the first part of the time span"""
    fo = (path / outfile).open("w")
    max_iter = 5
    t0 = 0.0
    t1 = 1.0
    rt = 0.0
    i = 0
    while i < max_iter:
        fo.write(f"Iteration {i}/{max_iter}: {t0:.2f} - {t1:.2f}\n")
        te, te0, te1 = detect_elastic_range(path, t0, t1, fo)
        rt = te / (t1 - t0)
        fo.write(
            f"{t0:.3f}  {t1:.3f}  {te:.3f} (e={int((te1-te0)/(t1-t0)*100)}%) {int(rt*100)}% \n"
        )
        if rt < 0.10:
            t1 /= 1.6
        elif rt > 0.45:
            t1 *= 2
        else:
            break
        i += 1


def detect_elastic_range(path, t0, t1, fo):
    """Performs several 1-timestep runs to find the boundary of the elastic range.
It receives the time range in which to look

path: base working path
t0: initial time
t1: final time
fo: file handler where to write iteration messages
"""

    param = json.loads((path / "ProjectParameters.json").resolve().read_text())

    #  Adapt params
    param["problem_data"]["echo_level"] = 0
    param["solver_settings"]["echo_level"] = 0
    param["solver_settings"]["line_search"] = False
    param["solver_settings"]["max_iteration"] = 1
    param["processes"]["my_processes"] = {}
    block = {
        "Parameters": {
            "model_part_name": "Microstructure.RVE",
            "filename": str(path / "is_elastic.tmp"),
        },
        "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
        "process_name": "IsInelastic",
        "python_module": "kratos_process_is_inelastic",
    }
    param["processes"]["my_processes"] = [block]
    #  Adapt model and material files
    #param["solver_settings"]["model_import_settings"]["input_filename"] = str(
    #    (path / "model").resolve()
    #)
    #param["solver_settings"]["material_import_settings"]["materials_filename"] = str(
    #    (path / "materials.json").resolve()
    #)

    # Begin iterations
    text = ""
    min_dt = (t1 - t0) / 10  # 10:~4 iters, 30:~5, 100:~7
    max_iter = 7
    i = 0
    while (t1 - t0) > min_dt and i < max_iter:
        t = t0 + (t1 - t0) / 2

        # Prepare case
        param["problem_data"]["end_time"] = t
        param["solver_settings"]["time_stepping"] = {}
        param["solver_settings"]["time_stepping"]["time_step"] = t
        # Run Kratos
        parameters = KratosMultiphysics.Parameters(json.dumps(param))
        model = KratosMultiphysics.Model()
        simulation = analysis.PBCAnalysis(model, parameters)
        simulation.Run()
        # Bound interval
        elastic = False
        if int((path / "is_elastic.tmp").resolve().read_text()) == 1:
            elastic = True

        if elastic:
            t0 = t
        else:
            t1 = t

        line = f"   {t0:.4f} - {t1:.4f} - Elastic: {elastic}\n"
        i += 1
        text += line
    fo.write(text)
    (path / "is_elastic.tmp").unlink()
    return t, t0, t1


def optimize_timesteps(path, t0, t1, te):

    #
    #       fix, 4 ts        fine, dt/2                 coarse, dt*2
    #   t0  .  .  .  . te  . ..t1..... . .   .   t2   .   .   .   .   tf
    #
    dt = (t1 - t0) / 40  # HARCODED, optimization parameter
    ts_table = [
        [t0, (te - t0) / 4],
        [te - (te - t0) * 0.2, (te - t0) / 4],
        [te + (t1 - te) * 0.0, dt / 2],
        [te + (t1 - te) * 0.4, dt * 2],
    ]

    param = json.loads((path / "ProjectParameters.json").resolve().read_text())
    param["problem_data"]["start_time"] = t0
    param["problem_data"]["end_time"] = t1
    param["solver_settings"]["time_stepping"] = {}
    param["solver_settings"]["time_stepping"]["automatic_time_step"] = False
    param["solver_settings"]["time_stepping"]["time_step_table"] = ts_table
    (path / "ProjectParameters.json").resolve().write_text(json.dumps(param, indent=4))


def print_histogram(values, bins=20, t0=0, t1=1):
    text = "Elastic range distribution:\n\n"
    B = [0] * bins
    for v in values:
        b = int(v / (t1 - t0) * bins)
        B[b] += 1
    for i, c in enumerate(B):
        text += "{:0.2f}: {:>3d} |{}\n".format(i * (t1 - t0) / bins, c, "*" * c)
    return text


def learn():
    cwd = Path.cwd()

    # Adjust elastic range by setting appropiate "tf"
    # if (case_path / "elastic.dat").exists() and not args["--force"]:
    adjust_elastic_range(cwd)

    # Adapt params
    line = (cwd / "elastic.dat").read_text().splitlines()[-1]
    t0, t1, te = [float(x) for x in line.split()[0:3]]
    optimize_timesteps(cwd, t0, t1, te)

    return


##############################
# MAIN
##############################

if __name__ == "__main__":
    learn()
