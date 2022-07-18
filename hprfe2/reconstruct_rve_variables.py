"""Reconstruct nodal and elementary fields.

Reconstructs total and fluctuant displacement (no damage nor stress),
which makes it compatible with micros and meso.

This is a workaround until we fix the reconstruction of ip fields in meso

Usage:
    reconstruct.py [-h]
    reconstruct.py [-v | -q] <root> <runtime_data>

Options:
-h --help     Show this
-v --verbose  Debug output
-q --quiet    Minimal output

Arguments:
root              Root path of the project
runtime_data      Generated run-time data file
"""
import os
from pathlib import Path
import logging
import math
import json
import numpy as np

# import h5py
from docopt import docopt

import meshio
from common import Common
import write_runtime_data
import KratosMultiphysics as KM

# import KratosMultiphysics.StructuralMechanicsApplication as SMA
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)

# import KratosMultiphysics.MultiscaleROMApplication

np.set_printoptions(
    linewidth=120,
    suppress=True,
)


##
## Functions for DAMAGE reconstruction
##


def q(r, E, yield_stress, inf_yield_stress, H0, H1):
    r0 = yield_stress / math.sqrt(E)
    q0 = r0  # strain_variable_init
    q1 = inf_yield_stress / math.sqrt(E)  # stress_variable_inf
    r1 = r0 + (q1 - q0) / H0
    if r < r0:
        return q0
    if r >= r0 and r < r1:
        return q0 + H0 * (r - r0)
    # Case r >= r1:
    return q1 + H1 * (r - r1)


def compute_elastic_tensor(E, NU):
    c_1 = E / ((1 + NU) * (1 - 2 * NU))
    c_2 = c_1 * (1 - NU)
    c_3 = c_1 * NU
    c_4 = c_1 * 0.5 * (1 - 2 * NU)
    elastic = np.zeros((6, 6))
    elastic[0, 0] = c_2
    elastic[0, 1] = c_3
    elastic[0, 2] = c_3
    elastic[1, 0] = c_3
    elastic[1, 1] = c_2
    elastic[1, 2] = c_3
    elastic[2, 0] = c_3
    elastic[2, 1] = c_3
    elastic[2, 2] = c_2
    elastic[3, 3] = c_4
    elastic[4, 4] = c_4
    elastic[5, 5] = c_4
    return elastic


def get_material_properties(model, props):
    material_properties = {}
    material_element_list = {}
    for m in props:
        material_name = m["model_part_name"]
        logger.debug("   - loading material {}".format(material_name))
        logger.debug(m["Material"]["Variables"])
        material_properties[material_name] = {}
        E = m["Material"]["Variables"]["YOUNG_MODULUS"]
        nu = m["Material"]["Variables"]["POISSON_RATIO"]
        yield_stress = m["Material"]["Variables"]["STRESS_LIMITS"][0]
        inf_yield_stress = m["Material"]["Variables"]["STRESS_LIMITS"][1]
        H0 = m["Material"]["Variables"]["HARDENING_PARAMETERS"][0]
        H1 = m["Material"]["Variables"]["HARDENING_PARAMETERS"][1]
        material_properties[material_name]["E"] = E
        material_properties[material_name]["nu"] = nu
        material_properties[material_name]["yield_stress"] = yield_stress
        material_properties[material_name]["inf_yield_stress"] = inf_yield_stress
        material_properties[material_name]["H0"] = H0
        material_properties[material_name]["H1"] = H1
        material_properties[material_name]["C"] = compute_elastic_tensor(E, nu)

        material_element_list[material_name] = []
        for elem in model[material_name].Elements:
            material_element_list[material_name].append(elem.Id)
    material_elem_map = {}
    for k, v in material_element_list.items():
        for idx in v:
            material_elem_map[idx] = k
    return material_properties, material_elem_map


##
## End functions for DAMAGE reconstruction
##


def strain_voigt_to_tensor(strain_vector):
    s_xx = strain_vector[0]
    s_yy = strain_vector[1]
    s_zz = strain_vector[2]
    s_xy = 0.5 * strain_vector[3]
    s_yz = 0.5 * strain_vector[4]
    s_xz = 0.5 * strain_vector[5]
    strain_tensor = np.array(
        [[s_xx, s_xy, s_xz], [s_xy, s_yy, s_yz], [s_yz, s_yz, s_zz]]
    )
    return strain_tensor


def get_mesh(rve_model):
    """Generate points and cells

    Returns:
        meshio.points -- points of the mesh
        meshio.cells -- cells of the mesh
    """
    mesh = meshio.read(rve_model)
    rve_cells = []
    for cell_block in mesh.cells:
        element_type = cell_block[0]
        # if "hexa" in element_type or "wedge" in element_type:
        if "line8" in element_type:
            rve_cells.append(meshio.CellBlock("hexahedron", cell_block[1]))
        if "line6" in element_type:
            rve_cells.append(meshio.CellBlock("wedge", cell_block[1]))
    return mesh.points, rve_cells


def init_kratos(mp_name, pmaterials, pmodel):
    """Load model and modelparts.

    Returns:
        dict -- for each element, location of beginning in the global dof vector
        dict -- for each element, number of integration points
    """
    parameters_dict = {
        "problem_data": {
            "problem_name": "High_Fidelity",
            "parallel_type": "OpenMP",
            "start_time": 0.0,
            "end_time": 0.99,
            "echo_level": 1,
        },
        "solver_settings": {
            "model_part_name": f"{mp_name}",
            "domain_size": 3,
            "echo_level": 1,
            "time_stepping": {},
            "solver_type": "Static",
            "model_import_settings": {
                "input_type": "mdpa",
                "input_filename": f"{pmodel}",
            },
            "material_import_settings": {"materials_filename": f"{pmaterials}"},
        },
    }

    model = KM.Model()
    parameters = KM.Parameters(json.dumps(parameters_dict))
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()
    return model, modelpart


def init_urt_data(model_part):
    # Select and initialize micros to write. AD-HOC. FIXME.
    eee = [
        (18, 2),
        (32, 0),
        (41, 0),
        (4, 0),
        (49, 0),
        (61, 3),
        (55, 3),
        (68, 5),
        (67, 3),
        (60, 5),
        (48, 7),
        (10, 3),
        (33, 0),
        (61, 7),
        (63, 7),
        (1, 2),
        (66, 4),
        (4, 7),
        (3, 1),
        (67, 5),
        (3, 3),
        (4, 5),
        (51, 6),
        (53, 7),
        (65, 2),
    ]
    uelems = [3, 0, 2, 4, 5]  # elements: 4, 18, x, 41, 49, 61
    uei_pairs = [eee[i] for i in uelems]

    urt_data = []
    for up in uei_pairs:
        ue = up[0]
        ui = up[1]
        filename = f"uruntime_{ue}_{ui}.json"
        self.urt_data.append(
            write_runtime_data.RuntimeData(filename, model_part, ue, ui)
        )
    return urt_data


class Reconstruct(Common):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        self.nr_voigt_comps = 6
        self.reconstruct_micro = False

    ##
    ## Functions for DAMAGE reconstruction
    ##

    def element_map(self):
        """Compute auxiliar vector with the index of an element in the global vector of dofs.

        Returns:
            dict -- for each element, location of beginning in the global dof vector
            dict -- for each element, number of integration points
        """
        count = 0
        elem_map = {}
        nips = {}
        for element in self.modelpart.Elements:
            elem_map[element.Id] = count
            nr_ip = len(
                element.CalculateOnIntegrationPoints(
                    KM.INTEGRATION_WEIGHT, self.modelpart.ProcessInfo
                )
            )
            nips[element.Id] = nr_ip
            count += nr_ip * self.nr_voigt_comps
        return elem_map, nips

    def reconstruc(self, runtime_data_path):
        # def initialize_runtime_micro(m, eips, nr_timesteps, nr_modes, nr_points):
        #    """Creates and initializes the micro runtime data file."""
        #    D = {}
        #    for x in m:
        #        eip = eips[x]
        #        e = eip[0]
        #        i = eip[1]
        #        d = {}
        #        d["multiscale_levels"] = 1
        #        d["nr_timesteps"] = nr_timesteps
        #        d["nr_modes"] = nr_modes
        #        d["nr_points"] = nr_points + 1
        #        d["interpolation_parameters"] = []
        #        d["macro_strain"] = []
        #        d["stress"] = []
        #        d["r_value"] = []
        #        D[x] = d
        #    return D

        # def append_runtime_micro(D, m, eips, strain_e, stress_e, c, rvalues):
        #    """Appends data to the micro runtime data file.
        #    eips is a list of tuples, containing the elemenet and the ip of the micro
        #    e.g. [(22, 0), (34, 7), (44, 3), ...]
        #    """
        #    for x in m:
        #        d = D[x]
        #        eip = eips[x]
        #        e = eip[0]
        #        i = eip[1]
        #        d["interpolation_parameters"].append(c[x])  # AD_HOC FIXME
        #        d["macro_strain"].append(list(strain_e[e][i, :]))
        #        d["stress"].append(list(stress_e[e][i, :]))
        #        d["r_value"].append(rvalues[x])
        #        D[x] = d
        #    return D

        # def write_runtime_micro(
        #    D,
        #    m,
        #    eips,
        # ):
        #    for x in m:
        #        e = eips[x][0]
        #        i = eips[x][1]
        #        write_json(f"runtime_{e}_{i}.json", D[x])

        # Load required data meso and micro
        logger.debug(f"Loading runtime data {runtime_data_path}")
        data = json.loads(runtime_data_path.read_text())
        nr_timesteps = data["nr_timesteps"]

        rve_interp_params = np.array(data[f"interpolation_parameters"])
        rve_macro_strain = np.array(data[f"macro_strain"])
        nr_modes = data[f"nr_modes"]
        nr_points = data[f"nr_points"] - 1
        logger.debug("Done")

        logger.debug("Loading strain bases")
        strain_modes = self.get_dataset("BASES", "STRAIN")[:, :nr_modes]
        logger.debug("Done")

        logger.debug("Loading strain correlation matrix")
        strain_correl = self.get_dataset("CORRELATION", "STRAIN", nr_modes)
        logger.debug("Done")

        logger.debug("Loading stress correlation matrix")
        stress_correl = self.get_dataset("CORRELATION", "STRESS", nr_modes, nr_points)
        logger.debug("Done")

        logger.debug("Loading rvalue correlation matrix")
        self.skip_damage_reconstruction = False
        try:
            r_value_correl = self.get_dataset(
                "CORRELATION", "RVALUE", nr_modes, nr_points
            )
            logger.debug("Done")
        except KeyError:
            logger.warning(
                "RVALUE correlation matrix not present. Skipping DAMAGE reconstruction"
            )
            self.skip_damage_reconstruction = True

        logger.debug("Loading rve data")
        rve_data = self.get_dataset("DATASET", "RVE", nr_modes, nr_points)

        logger.debug("Loading rve model")
        dset = self.get_dataset("TEMPLATE", "MODEL")
        logger.debug("Done")

        p_model = Path("model.mdpa")
        p_model.write_text(dset)
        rve_points, rve_cells = get_mesh(str(p_model))

        logger.debug("Loading rve materials")
        model_part_name = json.loads(
            self.get_dataset("TEMPLATE", "PARAMETERS_SAMPLING")
        )["solver_settings"]["model_part_name"]
        dset = self.get_dataset("TEMPLATE", "MATERIALS")
        p_materials = Path("materials.json")
        p_materials.write_text(dset)
        self.model, self.modelpart = init_kratos(
            model_part_name,
            str(p_materials.resolve()),
            str(p_model.resolve().parent / p_model.stem),
        )
        p_materials.unlink()
        p_model.unlink()

        # Get data from rve_data
        if not self.skip_damage_reconstruction:
            material_properties, material_elem_map = get_material_properties(
                self.model, rve_data["material_parameters"]["properties"]
            )

        ip_elem_map, nr_of_ips = self.element_map()

        if self.reconstruct_micro:
            init_urt_data()

        # Open XDMF file for writing field data for each timestep
        filename = "rve_reconstructed.xdmf"
        meshio.write_points_cells(filename, rve_points, rve_cells)
        with meshio.xdmf.TimeSeriesWriter(filename) as writer:
            writer.write_points_cells(rve_points, rve_cells)
            for t in range(nr_timesteps):
                logger.info("Timestep {}".format(t))

                logger.debug("Solving fluctuant displacement")
                displacement = np.dot(
                    strain_correl[:, :nr_modes], rve_interp_params[t, :]
                )
                displacement = np.reshape(displacement, (-1, 3))

                logger.debug("Solving total displacement")
                strain_macro = rve_macro_strain[t, :]
                strain_macro_tensor = strain_voigt_to_tensor(strain_macro)
                comp = np.dot(strain_macro_tensor, rve_points.T)
                total_displacement = comp.T + displacement

                stress_e, stress_h = self.compute_field_stress(
                    f"stress", stress_correl, data, t, nr_of_ips
                )
                strain_e, strain_h = self.compute_field_strain(
                    f"strain", strain_modes, rve_interp_params, t, nr_of_ips
                )

                ### Adding damage START
                if not self.skip_damage_reconstruction:
                    logger.debug("Solving damage")
                    damage_list = []
                    rvalue_list = []
                    r = np.dot(r_value_correl, data["r_value"][t])
                    r_in_elem = {}
                    for elem_id, nr_ips in nr_of_ips.items():
                        r_in_elem[elem_id] = r[:nr_ips]
                        r = r[nr_ips:]
                    for elem_id, nr_ips in nr_of_ips.items():
                        C = material_properties[material_elem_map[elem_id]]["C"]
                        E = material_properties[material_elem_map[elem_id]]["E"]
                        nu = material_properties[material_elem_map[elem_id]]["nu"]
                        yield_stress = material_properties[material_elem_map[elem_id]][
                            "yield_stress"
                        ]
                        inf_yield_stress = material_properties[
                            material_elem_map[elem_id]
                        ]["inf_yield_stress"]
                        H0 = material_properties[material_elem_map[elem_id]]["H0"]
                        H1 = material_properties[material_elem_map[elem_id]]["H1"]
                        r0 = yield_stress / math.sqrt(E)
                        damage = 0
                        rvalue = 0
                        for r in r_in_elem[elem_id]:
                            rvalue += r / nr_ips
                            if r < r0:
                                r = r0
                            d = 1 - q(r, E, yield_stress, inf_yield_stress, H0, H1) / r
                            damage += d / nr_ips
                        damage_list.append(damage)
                    element_damage = np.array(damage_list).reshape(
                        (-1, 1)
                    )  # formatting for meshio
                ### Adding damage END

                logger.debug("Writing timestep data")

                if self.reconstruct_micro:
                    uc = data["u_interpolation_parameters"][t]
                    ur = data["u_r_value"][t]
                    # mD = append_runtime_micro(
                    #    mD, micro_elems, eips, strain_e, stress_e, mc, mr
                    # )
                    append_urt_data(t, uc, strain_e, stress_e, ur)

                # Append XDMF Paraview data
                if self.skip_damage_reconstruction:
                    writer.write_data(
                        t,
                        point_data={
                            "DISPLACEMENT_FLUCT": np.reshape(displacement, (-1, 3)),
                            "DISPLACEMENT": total_displacement,
                        },
                        cell_data={
                            "STRAIN": strain_h,
                            "STRESS": stress_h,
                        },
                    )
                else:
                    writer.write_data(
                        t,
                        point_data={
                            "DISPLACEMENT_FLUCT": np.reshape(displacement, (-1, 3)),
                            "DISPLACEMENT": total_displacement,
                        },
                        cell_data={
                            "STRAIN": strain_h,
                            "STRESS": stress_h,
                            "DAMAGE": element_damage,
                        },
                    )

    def compute_field_stress(self, field, stress_correl, data, t, nr_of_ips):
        logger.debug(f"Computing {field} field")
        # Este es el stress en cada punto de gauss meso
        stress = np.dot(stress_correl, np.reshape(data[field][t], (-1, 1)))
        stress_r = stress.reshape((-1, 6))
        stress_e = {}
        # Este es el stress homogenizado en cada elemento meso
        stress_h = []
        idx = 0
        for elem_id, nr_ips in nr_of_ips.items():
            stress_e[elem_id] = stress_r[idx : idx + nr_ips, :]
            stress_h.append(np.mean(stress_e[elem_id], axis=0))
            idx += nr_ips
        return stress_e, stress_h

    def compute_field_strain(
        self, field, strain_modes, rve_interpolation_params, t, nr_of_ips
    ):
        logger.debug(f"Computing {field} field")
        strain_global = np.dot(strain_modes, rve_interpolation_params[t, :])
        strain_r = strain_global.reshape((-1, 6))
        strain_e = {}
        strain_h = []
        idx = 0
        for elem_id, nr_ips in nr_of_ips.items():
            strain_e[elem_id] = strain_r[idx : idx + nr_ips, :]
            strain_h.append(np.mean(strain_e[elem_id], axis=0))
            idx += nr_ips
        return strain_e, strain_h


#######################################
# Main
#######################################

# configure logger
# verbosity_level = logging.INFO
# if args.verbose:
#    verbosity_level = logging.DEBUG
verbosity_level = logging.DEBUG
logging.basicConfig(
    format="[%(asctime)s] %(message)s", datefmt="%H:%M:%S", level=verbosity_level
)
logger = logging.getLogger(__name__)

if __name__ == "__main__":

    ARGS = docopt(__doc__)

    RECONST = Reconstruct(root_path=Path(ARGS["<root>"]))
    RECONST.reconstruc(Path(ARGS["<runtime_data>"]))
