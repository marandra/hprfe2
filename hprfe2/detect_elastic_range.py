import sys
import json
from pathlib import Path
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis as periodic_bc_analysis
#import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis_xz as periodic_bc_analysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    # Este archivo leerlo de Common
    param= json.load(open("ProjectParameters.json"))

    # Adapt params: is_elastic process
    param["processes"]["my_processes"] = {}
    block = {
                "Parameters": {
                    "model_part_name": "Microstructure.RVE"
                },
                "kratos_module": "KratosMultiphysics.MultiscaleROMApplication",
                "process_name": "IsElastic",
                "python_module": "kratos_process_is_elastic"
            }
    param["processes"]["my_processes"] = [block]

    # Get params for learning
    t0 = param["problem_data"]["start_time"]
    t1 = param["problem_data"]["end_time"]

    # Compute params for binary search
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
        if int(Path("is_elastic.dat").read_text()) == 1:
            elastic = True

        if elastic:
            t0 = t
        else:
            t1 = t

        line = "{:.4f} - {:.4f} - Elastic: {}\n".format(t0, t1, elastic)
        i += 1
        text += line
    text += "{:.2f}\n".format(t)
    Path("elastic.dat").write_text(text)
    Path("is_elastic.dat").unlink()
