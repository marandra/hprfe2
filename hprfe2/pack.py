"""
PACK: pending description here.
"""

import json
from pathlib import Path
import numpy
import logging
from common import Common


logger = logging.getLogger(__name__)


def write_json(filename, data_dict):
    with open(filename, "w") as fo:
        json.dump(data_dict, fo, indent=2)


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def get_properties(rve_modelpart, iw_list):
    # read model materials
    out_prop = []
    for l in iw_list:
        elem_id = int(l[0])
        elem = rve_modelpart.GetElement(elem_id)
        prop_id = elem.Properties.Id
        out_prop.append(prop_id)
    return out_prop


def unpack_ip_data(iw_list):
    out_e = []
    out_ip = []
    out_w = []
    out_gip = []
    for l in iw_list:
        out_e.append(int(l[0]))
        out_ip.append(int(l[1]))
        out_w.append(float(l[2]))
        out_gip.append(int(l[3]))
    return out_e, out_ip, out_w, out_gip


def parse_strain_bases(strain_bases_filename, iw_list, nr_modes):
    strain_bases = numpy.load(str(strain_bases_filename), mmap_mode="r")
    strain_bases = strain_bases[:, :nr_modes]
    nr_comps = 6
    out_B = []
    for l in iw_list:
        gip = int(l[3])
        index = gip * nr_comps
        B = strain_bases[index : index + nr_comps, :]
        out_B.append(B.tolist())
    return out_B


def create_rve_params_structure(
    strain_bases_filename,
    rve_materials_filename,
    nr_modes,
    reduced_ip_set,
    rve_modelpart,
):
    """ gather and pack IP data for RVE constitutive law """
    rve_params = {}
    out_e, out_lip, out_w, out_gip = unpack_ip_data(reduced_ip_set)
    # required data
    rve_params["ip_global_id"] = out_gip
    rve_params["ip_weight"] = out_w
    rve_params["ip_property_id"] = get_properties(rve_modelpart, reduced_ip_set)
    rve_params["ip_strain_modes"] = parse_strain_bases(
        strain_bases_filename, reduced_ip_set, nr_modes
    )
    rve_params["material_parameters"] = read_json(rve_materials_filename)
    #  metadata
    rve_params["nr_modes"] = nr_modes
    rve_params["nr_reduced_ip"] = len(reduced_ip_set)
    rve_params["ip_element_id"] = out_e  # TODO: check if we use this data
    rve_params["ip_local_id"] = out_lip  # TODO: check if we use this data
    return rve_params


def write_datasets(common):
    """ docstring here """

    import KratosMultiphysics
    from KratosMultiphysics.StructuralMechanicsApplication import (
        structural_mechanics_analysis,
    )

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

    parameters = KratosMultiphysics.Parameters(json.dumps(parameters_dict))
    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters
    )
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()

    materials_fname = common.materials_fname
    for p in common.ip_subsets:
        roc_filename = common.bases_path / common.roc_fname(p)
        ip_set = numpy.loadtxt(roc_filename)
        for m in common.config["rve_data_modes"]:
            rve_fname = common.datasets_path / common.rve_fname(m, p)
            if common.skip_calculation(rve_fname):
                logger.info("File {} exists. Skipping calculation".format(rve_fname))
                continue
            logger.info("Generating {}".format(rve_fname))
            rve_params = create_rve_params_structure(
                common.get_bases_fname(common.config["strain_name"]),
                materials_fname,
                m,
                ip_set,
                modelpart,
            )
            write_json(rve_fname, rve_params)

    return


#######################################################################
#######################################################################

if __name__ == "__main__":

    import sys

    if len(sys.argv) > 1:
        C = Common(root_path=Path(sys.argv[1]))
    else:
        exit("Missing root_path argument.")

    write_datasets(C)
#    parameters_dict = {
#        "problem_data": {
#            "problem_name": "High_Fidelity",
#            "parallel_type": "OpenMP",
#            "start_time": 0.0,
#            "end_time": 0.99,
#            "echo_level": 1,
#        },
#        "solver_settings": {
#            "model_part_name": "Microstructure",
#            "domain_size": 3,
#            "echo_level": 1,
#            "time_stepping": {},
#            "solver_type": "Static",
#            "model_import_settings": {
#                "input_type": "mdpa",
#                "input_filename": "{}/model".format(co.training_path),
#            },
#            "material_import_settings": {
#                "materials_filename": "{}/materials.json".format(co.training_path)
#            },
#        },
#    }
#
#    parameters = KratosMultiphysics.Parameters(json.dumps(parameters_dict))
#    model = KratosMultiphysics.Model()
#    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
#        model, parameters
#    )
#    simulation.Initialize()
#    modelpart = simulation._GetSolver().GetComputingModelPart()
#
#    materials_fname = co.materials_fname
#    for p in co.ip_subsets:
#        roc_filename = co.roc_fname(p)
#        ip_set = numpy.loadtxt(roc_filename)
#        for m in co.context["rve_data_modes"]:
#            rve_fname = co.rve_fname(m, p)
#            if co.skip_calculation(rve_fname):
#                logger.info("File {} exists. Skipping calculation".format(rve_fname))
#                continue
#            logger.info("Generating {}".format(rve_fname))
#            rve_params = create_rve_params_structure(
#                co.get_bases_fname(co.context["strain_name"]),
#                materials_fname,
#                m,
#                ip_set,
#                modelpart,
#            )
#            write_json(rve_fname, rve_params)
