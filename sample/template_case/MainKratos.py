# makes KratosMultiphysics backward compatible with python 2.6 and 2.7
from __future__ import print_function, absolute_import, division

import sys
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication.periodic_bc_analysis as periodic_bc_analysis

"""
For user-scripting it is intended that a new class is derived
from StructuralMechanicsAnalysis to do modifications
"""

if __name__ == "__main__":

    fname = "ProjectParameters.json"
    if len(sys.argv) > 1:
        fname = sys.argv[1]

    with open(fname, "r") as parameter_file:
        parameters = KratosMultiphysics.Parameters(parameter_file.read())

    model = KratosMultiphysics.Model()
    simulation = periodic_bc_analysis.PBCAnalysis(model, parameters)
    simulation.Run()
