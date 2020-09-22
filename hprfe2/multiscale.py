"""
docstrings here
"""

import logging
import os
import json
from pathlib import Path
from common import Common


logger = logging.getLogger(__name__)


def create_properties_file(m_prop, c_prop, t_prop, quiet=False):
    """
    TODO: add docstrings here
    """
    test_props = json.loads(t_prop.read_text())
    strain_versor = test_props["processes"]["loads_process_list"][0]["Parameters"][
        "initial_strain"
    ]
    ampl = test_props["processes"]["loads_process_list"][0]["Parameters"][
        "lookuptable_mult"
    ][-1]

    model_props = json.loads(m_prop.read_text())
    # compute displacements u = E * x
    ss0, ss1, ss2, ss3, ss4, ss5 = strain_versor
    x0 = 1.0 * ss0 * ampl
    y0 = 0.5 * ss3 * ampl
    z0 = 0.5 * ss5 * ampl
    x1 = 0.5 * ss3 * ampl
    y1 = 1.0 * ss1 * ampl
    z1 = 0.5 * ss4 * ampl
    x2 = 0.5 * ss5 * ampl
    y2 = 0.5 * ss4 * ampl
    z2 = 1.0 * ss2 * ampl
    model_props["processes"]["list_boundary_processes"][1]["Parameters"]["value"] = [
        "{}*t".format(x0),
        "{}*t".format(x1),
        "{}*t".format(x2),
    ]
    model_props["processes"]["list_boundary_processes"][2]["Parameters"]["value"] = [
        "{}*t".format(y0),
        "{}*t".format(y1),
        "{}*t".format(y2),
    ]
    model_props["processes"]["list_boundary_processes"][3]["Parameters"]["value"] = [
        "{}*t".format(z0),
        "{}*t".format(z1),
        "{}*t".format(z2),
    ]

    if quiet:
        model_props["processes"]["my_processes"] = []
        model_props["output_processes"] = {}

    c_prop.write_text(json.dumps(model_props, indent=4))


def create_case_dir(rve, training, dataset):

    """
    Files and dirs structure:

    rve: root_path/multiscale/trajectory_35/_30m_400ip/
    source files: root_path/multiscale/MainKratos.py
                                       macro_model.mdpa
                                       ProjectParameters.json
    """

    # create dest dir
    rve.mkdir(parents=True, exist_ok=True)

    # adapt and copy materials file
    src = rve.parent.parent / "macro_materials.json"
    dest = rve / "macro_materials.json"
    rve_data_path = dataset / "rve{}.json".format(rve.name)
    materials = json.loads(src.read_text())
    materials["properties"][0]["Material"]["constitutive_law"]["Parameters"][
        "rve_data_filename"
    ] = str(rve_data_path.resolve())
    dest.write_text(json.dumps(materials, indent=4))

    # adapt and copy properties file
    m_prop = rve.parent.parent / "ProjectParameters.json"  # template properties file
    c_prop = rve / "ProjectParameters.json"  # destination case properties file
    t_prop = (
        training / rve.parent.name / "ProjectParameters.json"
    )  # test case properties file, get strain
    create_properties_file(m_prop, c_prop, t_prop)
    c_prop = rve / "ProjectParameters_quiet.json"
    create_properties_file(m_prop, c_prop, t_prop, quiet=True)

    # copy MainKratos.py
    src = rve.parent.parent / "MainKratos.py"
    dest = rve / "MainKratos.py"
    dest.write_text(src.read_text())

    # copy macro_model.mdpa
    src = rve.parent.parent / "macro_model.mdpa"
    dest = rve / "macro_model.mdpa"
    dest.write_text(src.read_text())

    return


def create_launch_script(case):
    """
    Writes temporary launch script for each case (run externally)
    """

    script_fname = "tmp_" + case.parent.name + case.name + ".bash"
    script = """\
export OMP_NUM_THREADS=1
export PYTHONPATH={}
export LD_LIBRARY_PATH={}
cd {}
/usr/bin/time -v -o time.dat python MainKratos.py ProjectParameters.json > outMainKratos
/usr/bin/time -v -o time_quiet.dat python MainKratos.py ProjectParameters_quiet.json > outMainKratos_quiet
cd {}
rm {}
""".format(
        os.environ["PYTHONPATH"],
        os.environ["LD_LIBRARY_PATH"],
        str(case),
        str(case.parent.parent),
        script_fname,
    )
    (case.parent.parent / script_fname).write_text(script)


def run(common):
    for c in common.config["validation_dataset"]:
        for m in common.config["rve_data_modes"]:
            for p in common.ip_subsets:
                rve_path = (
                    common.multiscale_path / common.case_name(c) / "_{}m_{}ip".format(m, p)
                ).resolve()
                create_case_dir(rve_path, common.training_path, common.datasets_path)
                create_launch_script(rve_path)
                logger.info("{} {}".format(rve_path.parent.name, rve_path.name))

#######################################
# main
#######################################

if __name__ == "__main__":

    import sys

    if len(sys.argv) > 1:
        common = Common(root_path=Path(sys.argv[1]))
    else:
        exit("Missing root_path argument.")

    run(common)
