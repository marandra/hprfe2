"""\

Usage:
    hprfe2 [-v] [-r PATH] validate

Options:
    -h --help             Print this help message and exit
    -v                    Verbose output (-v verbose, -vv debug)
    -r PATH --root=PATH   Specify the root path of the project, where the
                          configuration file must be located [default: .]

Writes and creates initial files structure for validation of computed bases.
"""

import logging
import os
import json
from pathlib import Path
from common import Common


logger = logging.getLogger(__name__)

MAIN = """
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)

with open("ProjectParameters.json", "r") as fp:
    parameters = KratosMultiphysics.Parameters(fp.read())
model = KratosMultiphysics.Model()
simulation = StructuralMechanicsAnalysis(model, parameters)
simulation.Run()
"""

MATERIAL = """
{  
    "properties": [{
        "model_part_name": "Structure.MATERIAL_MULTISCALE",
        "properties_id": 1,
        "Material": {
            "name": "multiscale",
            "constitutive_law": {
                "name": "RVELaw",
                "Parameters": {
                        "rve_data_filename": "to_be_filled_by_script",
                        "convergence_criterion": "residual_criterion",
                        "residual_relative_tolerance": 1e-2,
                        "residual_absolute_tolerance": 0.0,
                        "max_iteration": 20,
                        "verbose": 1
                        }
                },
            "Variables": {},
            "Tables":  {}
            }
       }]
}
"""

MODEL = """
Begin ModelPartData
End ModelPartData

Begin Properties 0
End Properties
Begin Nodes
    1		  0.0	   0.0	   0.0
    2		  1.0	   0.0	   0.0
    3  	  0.0	   1.0	   0.0
    4		  0.0	   0.0	   1.0
End Nodes

Begin Elements SmallDisplacementElement3D4N
    1          0         4	 3	 2	 1 
End Elements

Begin SubModelPart MATERIAL_MULTISCALE
    Begin SubModelPartNodes
    End SubModelPartNodes
    Begin SubModelPartElements
    1
    End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart

Begin SubModelPart PINNED
    Begin SubModelPartNodes
        1
    End SubModelPartNodes
    Begin SubModelPartElements
    End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart

Begin SubModelPart DISPL_X
    Begin SubModelPartNodes
        2
    End SubModelPartNodes
    Begin SubModelPartElements
    End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart

Begin SubModelPart DISPL_Y
    Begin SubModelPartNodes
        3
    End SubModelPartNodes
    Begin SubModelPartElements
    End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart

Begin SubModelPart DISPL_Z
    Begin SubModelPartNodes
        4
    End SubModelPartNodes
    Begin SubModelPartElements
    End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart

Begin SubModelPart MACRO
    Begin SubModelPartNodes
         1
         2
         3
         4
    End SubModelPartNodes
    Begin SubModelPartElements
         1
    End SubModelPartElements
    Begin SubModelPartConditions
    End SubModelPartConditions
End SubModelPart
"""

PARAMS = """
{
    "problem_data": {
        "problem_name": "Multiscale",
        "parallel_type": "OpenMP",
        "start_time": 0.0,
        "end_time": 0.99,
        "echo_level": 1
    },
    "solver_settings": {
        "model_part_name": "Structure",
        "domain_size": 3,
        "echo_level": 1,
        "time_stepping": {
            "time_step": 0.025
        },
        "solver_type": "Static",
        "model_import_settings": {
            "input_type": "mdpa",
            "input_filename": "macro_model"
        },
        "material_import_settings": {
            "materials_filename": "macro_materials.json"
        },
        "line_search": true,
        "use_old_stiffness_in_first_iteration": true,
        "move_mesh_flag": false,
        "convergence_criterion": "residual_criterion",
        "residual_relative_tolerance": 1e-2,
        "residual_absolute_tolerance": 0.0,
        "max_iteration": 10,
        "linear_solver_settings": {
            "solver_type": "amgcl",
            "krylov_type": "cg",
            "max_iteration": 500,
            "scaling": false,
            "verbosity": 1
        },
        "rotation_dofs": false,
        "compute_reactions": true,
        "auxiliary_variables_list": []
        },
    "processes": {
        "my_processes": [{
           "python_module": "kratos_process_write_rve_reconstruction_data",
           "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
           "process_name": "WriteRveReconstructionData",
           "Parameters": {
               "model_part_name": "Structure.MACRO",
               "filename": "rve_runtime_data_el1_ip0.json",
               "element": 1,
               "integration_point": 0
               }
           },{
            "python_module": "kratos_process_write_elements_homogenized_output",
            "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
            "process_name": "WriteElementsHomogenizedOutput",
            "Parameters": {
                "model_part_name": "Structure.MACRO",
                "filename": "homogenized_stress.dat",
                "variable_name": "CAUCHY_STRESS_VECTOR"
                }
            } 
        ],
        "list_initial_processes": [],
        "list_boundary_processes": [{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {
                "model_part_name": "Structure.PINNED",
                "variable_name": "DISPLACEMENT",
                "constrained": [true, true, true],
                "value": [0.0, 0.0, 0.0],
                "interval": [0.0, "End"]
                }
        },{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {
                "model_part_name": "Structure.DISPL_X",
                "variable_name": "DISPLACEMENT",
                "constrained": [true, true, true],
                "value": ["0.0*t", "0.0*t", "0.0*t"],
                "interval": [0.0, "End"]
                }
            },{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {
                "model_part_name": "Structure.DISPL_Y",
                "variable_name": "DISPLACEMENT",
                "constrained": [true, true, true],
                "value": ["0.0*t", "0.0*t", "0.0*t"],
                "interval": [0.0, "End"]
                }
            },{
            "python_module": "assign_vector_variable_process",
            "kratos_module": "KratosMultiphysics",
            "process_name": "AssignVectorVariableProcess",
            "Parameters": {
                "model_part_name": "Structure.DISPL_Z",
                "variable_name": "DISPLACEMENT",
                "constrained": [true, true, true],
                "value": ["0.0*t", "0.0*t", "0.0*t"],
                "interval": [0.0, "End"]
                }
            }],
        "loads_process_list": []
    },
    "output_processes" : {},
    "restart_options": {
        "SaveRestart": false,
        "RestartFrequency": 0,
        "LoadRestart": false,
        "Restart_Step": 0
    },
    "constraints_data": {
        "incremental_load": false,
        "incremental_displacement": false
    }
}
"""


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

    # Create base directory
    p = common.multiscale_path
    if not p.exists():
        p.mkdir()
        logger.info("Created directory {}".format(p))

    # Write template files to validation directory
    dest = common.multiscale_path / "macro_model.mdpa"
    dest.write_text(MODEL)
    dest = common.multiscale_path / "macro_materials.json"
    dest.write_text(MATERIAL)
    dest = common.multiscale_path / "ProjectParameters.json"
    dest.write_text(PARAMS)
    dest = common.multiscale_path / "MainKratos.py"
    dest.write_text(MAIN)
    logger.info("Written template files")

    # Create case structure
    for c in common.config["validation_dataset"]:
        for m in common.config["rve_data_modes"]:
            for p in common.ip_subsets:
                rve_path = (
                    common.multiscale_path / common.case_name(c) / "_{}m_{}ip".format(m, p)
                ).resolve()
                create_case_dir(rve_path, common.training_path, common.datasets_path)
                create_launch_script(rve_path)
                logger.info("{} {}".format(rve_path.parent.name, rve_path.name))
