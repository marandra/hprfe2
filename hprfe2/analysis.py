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

import json
from pathlib import Path
import sys

import KratosMultiphysics
import KratosMultiphysics.MultiscaleROMApplication
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication import (
    structural_mechanics_analysis,
)

RVE0_CONST = 0.030e-3  # 1 iteration, evaluating 1 RVE0
RVE1_CONST = 0.020e-6  #  const for solving system
NL_ITER = 6  # average non-linear iterations and substep line search


class Material:
    """doc"""

    def __init__(self, value, count, nmodes=0, children=None):
        self.prop = value
        self.count = count
        self.nmodes = nmodes
        self.children = children

    def __repr__(self, level=0):
        # tabulation
        # ret = f"{self.estimate_time():0.6f}s"
        ret = ""
        ret += "    " * level

        # quantity, id, part name
        pid = self.prop["properties_id"]
        ret += f"{self.count:6d}x id {pid:3d}: {self.prop['model_part_name']} "

        # nr modes and nr points
        if self.nmodes != 0:
            nip = 0
            for child in self.children:
                nip += child.count
            ret += f"({self.nmodes} modes, {nip} points)"

        # eol
        ret += "\n"

        # recurse into children
        if self.children is not None:
            for child in self.children:
                ret += child.__repr__(level + 1)

        return ret

    def estimate_time(self):
        """docs"""
        children_time = 0.0
        for child in self.children:
            children_time += child.estimate_time()
        unit_time = RVE0_CONST + NL_ITER * (
            RVE1_CONST * self.nmodes ** 2 + children_time
        )
        return self.count * unit_time


def analyze(props, count):
    """Add docstring"""
    materials = []
    for prop in props:
        nmodes = 0
        children = []
        name = prop["Material"]["constitutive_law"]["name"]
        if "RVELaw" in name:
            rve_fname = prop["Material"]["constitutive_law"]["Parameters"][
                "rve_data_filename"
            ]
            rve_props, rve_count, nmodes = get_properties_from_rve(rve_fname)
            children = analyze(rve_props, rve_count)

        materials.append(Material(prop, count[prop["properties_id"]], nmodes, children))
        # else:
        #    for k, v in prop['Material']['Variables'].items():
        #        print(f"{offset}   {k}: {v}")
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
    nmodes = len(rve["ip_strain_modes"][0][0])
    return props, count, nmodes


def load_case(case):
    """displacement docstrings here"""
    params = json.loads((case / "ProjectParameters.json").read_text())

    # remove processes
    params["processes"]["my_processes"] = []
    # make paths absolute
    model_p = Path(params["solver_settings"]["model_import_settings"]["input_filename"])
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
        nip = len(elem.GetIntegrationPoints())
        idx.extend([elem.Properties.Id] * nip)
    for i in set(idx):
        count[i] = idx.count(i)
    # simulation.RunSolutionLoop()
    # simulation.Finalize()
    props = json.loads(materials_p.read_text())["properties"]
    return props, count


def run(case):
    """Add docstring"""
    properties, count = load_case(case)
    materials = analyze(properties, count)
    print()
    print("-------------------------------------------------------------")
    print()
    print("Materials structure")
    for material in materials:
        print(material)
    time = 0
    for material in materials:
        time += material.estimate_time()
    print("Estimated times")
    print(f" - one iteration: {time:0.6f}s")
    print(
        f" - trajectory (40 steps), no output:  {40 * time:0.6f}s - {2 * 40 * time:0.6f}s"
    )
    print(f" - trajectory with output (validation): {1.3 * 40 * time:0.6f}s")
    print(
        f" - trajectory with output (sampling): {2 * 40 * time:0.6f}s - {3 * 40 * time:0.6f}s"
    )


####################3

if __name__ == "__main__":
    run(Path(sys.argv[1]))
