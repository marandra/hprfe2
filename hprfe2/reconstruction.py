"""
RECONSTRUCTION: description here
"""

import logging
import json
import numpy
import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import (
    MechanicalSolver,
)
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)

logger = logging.getLogger(__name__)

#
# Classes for setting Kratos solvers
#
class DisplacementReconstructionSolver(MechanicalSolver):
    def __init__(self, model, custom_settings):
        super(DisplacementReconstructionSolver, self).__init__(model, custom_settings)
        KratosMultiphysics.Logger.PrintInfo(
            "::[Displacement Reconstruction Solver]:: ", "Construction finished"
        )

    def _create_solution_scheme(self):
        return KratosMultiphysics.ResidualBasedIncrementalUpdateStaticScheme()

    def _create_builder_and_solver(self):
        linear_solver = self.get_linear_solver()
        builder_and_solver = KratosMultiphysics.MultiscaleROMApplication.ResidualBasedBlockBuilderAndSolverCustom(
            linear_solver
        )
        return builder_and_solver


class DisplacementReconstructionAnalysis(AnalysisStage):
    def __init__(self, model, project_parameters):
        super(DisplacementReconstructionAnalysis, self).__init__(
            model, project_parameters
        )

    #### Must be defined ####
    def _CreateSolver(self):
        solver = DisplacementReconstructionSolver(
            self.model, self.project_parameters["solver_settings"]
        )
        return solver


#
# Functions for DISPLACEMENT reconstruction
#
def reconstruct_displacement_all(common):
    ### workaround until we can make Kratos change the elemen type
    global_index_path = common.bases_path / "auxiliar_global_index"
    model_custom_path = common.bases_path / "model_custom.mdpa"
    model_original_text = (common.training_path / "model.mdpa").read_text()
    model_custom_text = model_original_text.replace(
        "DisplacementElement", "DisplacementCustomElement"
    )
    model_custom_path.write_text(model_custom_text)

    for pair in common.config["reconstruction_pairs"]:
        nm = pair[0]
        if common.has_dataset("CORRELATION", "STRAIN", nm):
            logger.info(
                f'CORRELATION {common.name_dataset("STRAIN", nm)} exists. Skipping.'
            )
            continue
        else:
            reconstruct_displacement(common, nm)

    global_index_path.unlink(missing_ok=True)
    return


def reconstruct_displacement(common, n_modes):
    """displacement docstrings here"""

    # Define parameters for reconstruction
    model_custom_path = common.bases_path / "model_custom.mdpa"
    global_index_path = common.bases_path / "auxiliar_global_index"
    materials_fname = common.training_path / "materials.json"

    params_reconstr_dict = {
        "problem_data": {
            "problem_name": "High_Fidelity",
            "parallel_type": "OpenMP",
            "start_time": 0.0,
            "end_time": 0.99,
            "echo_level": 1,
        },
        "solver_settings": {
            "model_part_name": "Microstructure",
            "domain_size": 3,
            "echo_level": 1,
            "time_stepping": {"time_step": 1.0},
            "solver_type": "Static",
            "model_import_settings": {
                "input_filename": str(model_custom_path.with_suffix("")),
                "input_type": "mdpa",
            },
            "material_import_settings": {
                "materials_filename": str(materials_fname),
            },
            "linear_solver_settings": {
                "solver_type": "amgcl",
                "krylov_type": "cg",
                "max_iteration": 500,
                "scaling": False,
                "verbosity": 1,
            },
            "line_search": False,
            "convergence_criterion": "residual_criterion",
            "residual_relative_tolerance": 1e-4,
            "residual_absolute_tolerance": 0.0,
            "max_iteration": 1,
            "rotation_dofs": False,
            "compute_reactions": False,
            "move_mesh_flag": False,
            "block_builder": True,
            "auxiliary_variables_list": [],
        },
        "processes": {
            "my_processes": [
                {
                    "python_module": "kratos_process_load_modes_to_properties",
                    "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
                    "process_name": "LoadModesToProperties",
                    "Parameters": {
                        "model_part_name": "Microstructure.RVE",
                        "global_index_filename": str(global_index_path),
                        "number_modes_to_load": n_modes,
                        "root_path": str(common.root_path.resolve()),
                    },
                }
            ],
            "list_initial_processes": [],
            "list_boundary_processes": [
                {
                    "python_module": "assign_vector_variable_process",
                    "kratos_module": "KratosMultiphysics",
                    "process_name": "AssignVectorVariableProcess",
                    "Parameters": {
                        "model_part_name": "Microstructure.PINNED",
                        "variable_name": "DISPLACEMENT",
                        "constrained": [True, True, True],
                        "value": [0.0, 0.0, 0.0],
                        "interval": [0.0, "End"],
                    },
                }
            ],
            "loads_process_list": [],
        },
        "output_processes": {},
        "restart_options": {
            "SaveRestart": False,
            "RestartFrequency": 0,
            "LoadRestart": False,
            "Restart_Step": 0,
        },
        "constraints_data": {
            "incremental_load": False,
            "incremental_displacement": False,
        },
    }

    parameters_reconstr = KratosMultiphysics.Parameters(
        json.dumps(params_reconstr_dict)
    )

    #  Generate auxiliar data structure
    parameters_dict = {
        "problem_data": {
            "problem_name": "High_Fidelity",
            "parallel_type": "OpenMP",
            "start_time": 0.0,
            "end_time": 0.99,
            "echo_level": 1,
        },
        "solver_settings": {
            "model_part_name": "Microstructure",
            "domain_size": 3,
            "echo_level": 1,
            "time_stepping": {},
            "solver_type": "Static",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": "{}/model".format(common.training_path),
            },
            "material_import_settings": {
                "materials_filename": "{}/materials.json".format(common.training_path)
            },
        },
    }
    parameters_aux = KratosMultiphysics.Parameters(json.dumps(parameters_dict))
    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters_aux
    )
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()
    for elem in modelpart.Elements:
        nr_comp = len(
            elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.STRAIN, modelpart.ProcessInfo
            )[0]
        )
        break
    idx_vector = []
    count = 0
    for elem in modelpart.Elements:
        idx_vector.append(count)
        nr_ips = len(
            elem.CalculateOnIntegrationPoints(
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, modelpart.ProcessInfo
            )
        )
        count = count + nr_ips * nr_comp
    fname = parameters_reconstr["processes"]["my_processes"][0]["Parameters"][
        "global_index_filename"
    ].GetString()
    with open(fname, "w") as ofile:
        for idx in idx_vector:
            ofile.write("{}\n".format(idx))
    # end of generating auxiliar file

    # Reconstruction
    model = KratosMultiphysics.Model()
    simulation = DisplacementReconstructionAnalysis(model, parameters_reconstr)
    # we replace .Run() by the code below so we can remove conditions
    # (and in the future replace elements, no we don need to modify model.mdpa)
    # simulation.Run()
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()

    for condition in modelpart.Conditions:
        condition.Set(KratosMultiphysics.TO_ERASE)
    modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)

    # settings = KratosMultiphysics.Parameters("""
    #    {
    #        "element_name": "SmallDisplacementCustomElement3D8N",
    #        "condition_name": ""
    #    }
    #    """)
    # KratosMultiphysics.ReplaceElementsAndConditionsProcess(modelpart, settings).Execute()
    simulation.RunSolutionLoop()
    simulation.Finalize()

    global_index_path.unlink(missing_ok=True)

    return


#
# Functions for DAMAGE Reconstruction
#
def load_rve_data(rve_data):
    logger.debug("Reading reduced set integration points")
    reduced_ip_set = rve_data["ip_global_id"]
    logger.debug("Nr ip detected: {}".format(numpy.shape(reduced_ip_set)[0]))
    logger.debug("Reading reduced set integration weights")
    reduced_ip_weights = numpy.array(rve_data["ip_weight"])
    logger.debug("Nr weights detected: {}".format(numpy.shape(reduced_ip_weights)[0]))
    return reduced_ip_set, reduced_ip_weights


def load_energy_modes(common, reduced_ip_set, nr_modes):
    modes = common.get_dataset("BASES", "RVALUE")[:, :nr_modes]
    reduced_modes = modes[reduced_ip_set, :]
    logger.info(
        "Modes matrix {} {} - Reduced modes matrix: {} {}".format(
            numpy.shape(modes)[0],
            numpy.shape(modes)[1],
            numpy.shape(reduced_modes)[0],
            numpy.shape(reduced_modes)[1],
        )
    )
    return modes, reduced_modes


def compute_reconstruction_system(
    reduced_ip_weights, energy_modes, reduced_energy_modes
):
    logger.info("Computing COMPLETE system")
    logger.debug("-- A = reduced modes.T * weights * reduced modes")
    reduced_ip_weights_diag = numpy.diag(reduced_ip_weights)

    weighted_reduced_modes_transposed = numpy.dot(
        reduced_energy_modes.T, reduced_ip_weights_diag
    )
    A = numpy.dot(weighted_reduced_modes_transposed, reduced_energy_modes)

    logger.debug("-- checking A is not singular")
    rankA = numpy.linalg.matrix_rank(A)
    logger.debug("A: {}".format(numpy.shape(A)))
    logger.debug("rank A: {}".format(numpy.linalg.matrix_rank(A)))
    if rankA != numpy.shape(A)[0]:
        logger.info("Matrix rank not complete (Too many ROC points?). Aborting.")
        exit()
    logger.debug("-- inverse A")
    Ainv = numpy.linalg.inv(A)

    logger.debug("-- modes * invA * modes.T * weights ")
    aux_1 = numpy.dot(Ainv, weighted_reduced_modes_transposed)
    aux_2 = numpy.dot(energy_modes, aux_1)
    return aux_2


def compute_system(common, rve_data, nr_modes):
    """docstrings"""
    reduced_ip_set, reduced_ip_weights = load_rve_data(rve_data)
    modes, reduced_modes = load_energy_modes(common, reduced_ip_set, nr_modes)
    A = compute_reconstruction_system(reduced_ip_weights, modes, reduced_modes)
    return A


def reconstruct_damage_all(common):
    """Computes data necessary for later reconstruction of the damage.
    Skips computation if file exists and option 'reuse_existing_file' is set."""
    for pair in common.config["reconstruction_pairs"]:
        nm = pair[0]
        np = pair[1]
        if common.has_dataset("CORRELATION", "RVALUE", nm, np):
            logger.info(
                f'CORRELATION {common.name_dataset("RVALUE", nm, np)} exists. Skipping.'
            )
            continue
        else:
            reconstruct_damage(common, nm, np)
    return


def reconstruct_damage(common, nm, np):
    rve_data = common.get_dataset("DATASET", "RVE", nm, np)
    logger.info(f'Computing CORRELATION {common.name_dataset("RVALUE", nm, np)}')
    A = compute_system(common, rve_data, nm)
    common.set_dataset(A, "CORRELATION", "RVALUE", nm, np)
    return
