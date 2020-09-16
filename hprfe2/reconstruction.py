"""
RECONSTRUCTION: description here
"""

from pathlib import Path
import json

import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.analysis_stage import AnalysisStage
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_solver import (
    MechanicalSolver,
)
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)
from common import Common


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


def reconstruct_displacement(common):
    """ docstrings here """

    # Define parameters for reconstruction
    model_custom_fname = common.bases_path / "model_custom.mdpa"
    strain_bases_fname = common.get_bases_fname("STRAIN")
    n_modes = 30

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
                "input_filename": str(model_custom_fname.with_suffix("")),
                "input_type": "mdpa",
            },
            "material_import_settings": {
                "materials_filename": str(common.training_path / "materials.json"),
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
                        "modes_filename": str(strain_bases_fname),
                        "modes_file_format": "binary",
                        "global_index_filename": str(
                            common.bases_path / "auxiliar_global_index"
                        ),
                        "number_modes_to_load": n_modes,
                        "modes_to_nodes_matrix_filename": str(
                            common.bases_path
                            / "correlation_strain_{}.npy".format(n_modes)
                        ),
                        "modes_to_nodes_matrix_file_format": "binary",
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

    # workaround until we can make Kratos change the elemen type
    model_original_text = (common.training_path / "model.mdpa").read_text()
    model_custom_text = model_original_text.replace(
        "DisplacementElement", "DisplacementCustomElement"
    )
    model_custom_fname.write_text(model_custom_text)

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
                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, modelpart.ProcessInfo
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

    return


###############################################################
###############################################################

if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1:
        co = Common(root_path=Path(sys.argv[1]))
    else:
        exit("Missing root_path argument.")
#
#    # Read parametres for reconstruction
#    with open("../ProjectParameters_correlation.json", "r") as parameter_file:
#        parameters_reconstr = KratosMultiphysics.Parameters(parameter_file.read())
#
#    #  Generate auxiliar data structure
#    parameters_dict = {
#        "problem_data": {
#        "problem_name": "High_Fidelity",
#        "parallel_type": "OpenMP",
#        "start_time": 0.0,
#        "end_time": 0.99,
#        "echo_level": 1,
#    },
#    "solver_settings": {
#        "model_part_name": "Microstructure",
#        "domain_size": 3,
#        "echo_level": 1,
#        "time_stepping": {},
#        "solver_type": "Static",
#        "model_import_settings": {
#            "input_type": "mdpa",
#            "input_filename": "{}/model".format(co.training_path),
#        },
#        "material_import_settings": {
#            "materials_filename": "{}/materials.json".format(co.training_path)
#        },
#    },
#    }
#    parameters_aux = KratosMultiphysics.Parameters(json.dumps(parameters_dict))
#    model = KratosMultiphysics.Model()
#    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
#        model, parameters_aux
#    )
#    simulation.Initialize()
#    modelpart = simulation._GetSolver().GetComputingModelPart()
#    for elem in modelpart.Elements:
#        nr_comp = len(
#            elem.CalculateOnIntegrationPoints(
#                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, modelpart.ProcessInfo
#            )[0]
#        )
#        break
#    idx_vector = []
#    count = 0
#    for elem in modelpart.Elements:
#        idx_vector.append(count)
#        nr_ips = len(
#            elem.CalculateOnIntegrationPoints(
#                KratosMultiphysics.GREEN_LAGRANGE_STRAIN_VECTOR, modelpart.ProcessInfo
#            )
#        )
#        count = count + nr_ips * nr_comp
#    fname = parameters_reconstr["processes"]["my_processes"][0]["Parameters"][
#        "global_index_filename"
#    ].GetString()
#    with open(fname, "w") as ofile:
#        for idx in idx_vector:
#            ofile.write("{}\n".format(idx))
#    # end of generating auxiliar file
#
#    # Reconstruction
#    model = KratosMultiphysics.Model()
#    simulation = DisplacementReconstructionAnalysis(model, parameters_reconstr)
#    # we replace .Run() by the code below so we can remove conditions
#    # (and in the future replace elements, no we don need to modify model.mdpa)
#    # simulation.Run()
#    simulation.Initialize()
#    modelpart = simulation._GetSolver().GetComputingModelPart()
#
#    for condition in modelpart.Conditions:
#        condition.Set(KratosMultiphysics.TO_ERASE)
#    modelpart.RemoveConditionsFromAllLevels(KratosMultiphysics.TO_ERASE)
#
#    # settings = KratosMultiphysics.Parameters("""
#    #    {
#    #        "element_name": "SmallDisplacementCustomElement3D8N",
#    #        "condition_name": ""
#    #    }
#    #    """)
#    # KratosMultiphysics.ReplaceElementsAndConditionsProcess(modelpart, settings).Execute()
#
#    simulation.RunSolutionLoop()
#    simulation.Finalize()
