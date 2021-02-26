"""\

Usage:
    hprfe2 [-v] [-r PATH] summary -c case_path

Options:
    -h --help             Print this help message and exit
    -v                    Verbose output (-v verbose, -vv debug)
    -r PATH --root=PATH   Specify the root path of the project, where the
                          configuration file must be located [default: .]
    -c --case CASE_PATH   Path of case to analyze

Summary of the model
"""

import logging
import json
from pathlib import Path
import sys
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
import KratosMultiphysics.MultiscaleROMApplication
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)

logger = logging.getLogger(__name__)


class Material():
    """doc"""
    def __init__(self, value, children):
        self.prop = value
        self.count = 0
        self.children = children

    def __repr__(self, level=0):
        #line = f"level {level}{offset} "
        #line += f"{count[pid]:6d}x "
        #line += f"property {pid:3d} {rve_nmodes} '{prop['Material']['name']}'"
        #line += f" - {clname}"
        ret = ""
        pid = self.prop["properties_id"]
        #ret += "\t" * level + f"{self.count[pid]:6d}x id {pid:3d} {prop['model_part_name']}" + "\n"
        ret += "\t" * level + f"{self.count:6d}x id {pid:3d} {self.prop['model_part_name']}" + "\n"
        if self.children is not None:
            for child in self.children:
                ret += child.__repr__(level+1)
        return ret

def analyze(props):
    """Add docstring"""
    #rve_nmodes = -1
    #print(material.count)
    #print(material.props)
    materials = []
    for prop in props:
        children = []
        name = prop["Material"]["constitutive_law"]["name"]
        if "RVELaw" in name:
            rve_fname = prop["Material"]["constitutive_law"]["Parameters"][ "rve_data_filename" ]
            rve_props, rve_count, rve_nmodes = get_properties_from_rve(rve_fname)
            children = analyze(rve_props)

        materials.append(Material(prop, children))
        # else:
        #    for k, v in prop['Material']['Variables'].items():
        #        print(f"{offset}   {k}: {v}")
        #    print()
    return materials


def get_properties_from_rve(rve_fname):
    """Add docstring"""
    rve = json.load(open(rve_fname, "r"))
    props = rve["material_parameters"]["properties"]
    count = {}
    for i in set(rve["ip_property_id"]):
        count[i] = 0
    for i in rve["ip_property_id"]:
        count[i] += 1
    return props, count, 10


def summary_material(props, count, level):
    """Add docstring"""
    level -= 1
    offset = "    " * level
    rve_nmodes = -1
    for prop in props:
        clname = prop["Material"]["constitutive_law"]["name"]
        if "RVELaw" in clname:
            rve_fname = prop["Material"]["constitutive_law"]["Parameters"][ "rve_data_filename" ]
            rve_props, rve_count, rve_nmodes = get_properties_from_rve(rve_fname)
            level = summary_material(rve_props, rve_count, level)
        # else:
        #    for k, v in prop['Material']['Variables'].items():
        #        print(f"{offset}   {k}: {v}")
        #    print()
        pid = prop["properties_id"]
        line = f"level {level}{offset} "
        if count is not None:
            line += f"{count[pid]:6d}x "
        line += f"property {pid:3d} {rve_nmodes} '{prop['Material']['name']}'"
        line += f" - {clname}"
        print(line)
    return level + 1


def load_case(case):
    """displacement docstrings here """
    params = json.loads((case / "ProjectParameters.json").read_text())

    # make paths absolute
    model_p = Path( params["solver_settings"]["model_import_settings"]["input_filename"])
    if not model_p.is_absolute():
        model_p = case / model_p
    params["solver_settings"]["model_import_settings"]["input_filename"] = str(model_p)
    materials_p = Path(
        params["solver_settings"]["material_import_settings"]["materials_filename"]
    )
    if not materials_p.is_absolute():
        materials_p = case / materials_p
    params["solver_settings"]["material_import_settings"]["materials_filename"] = str(
        materials_p
    )

    parameters_aux = KratosMultiphysics.Parameters(json.dumps(params))
    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters_aux
    )
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()
    # cl_descr = modelpart.Properties[1].GetValue(KratosMultiphysics.CONSTITUTIVE_LAW)
    idx = []
    count = {}
    for elem in modelpart.Elements:
        idx.append(elem.Properties.Id)
    for i in set(idx):
        count[i] = idx.count(i)
    # simulation.RunSolutionLoop()
    # simulation.Finalize()
    props = json.loads(materials_p.read_text())["properties"]
    return props, count


def run(case):
    """Add docstring"""
    properties, count = load_case(case)
    materials = analyze(properties)
    print()
    print("Materials structure:")
    #summary_material(properties, count, 3)
    for material in materials:
        print(material)



####################3

if __name__ == "__main__":
    run(Path(sys.argv[1]))
