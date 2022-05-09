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
import numpy
import h5py
from docopt import docopt
import meshio
from common import Common
import KratosMultiphysics
import KratosMultiphysics.StructuralMechanicsApplication
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)
import KratosMultiphysics.MultiscaleROMApplication

numpy.set_printoptions(
    linewidth=120,
    suppress=True,
)


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
    elastic = numpy.zeros((6, 6))
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


def strain_voigt_to_tensor(strain_vector):
    s_xx = strain_vector[0]
    s_yy = strain_vector[1]
    s_zz = strain_vector[2]
    s_xy = 0.5 * strain_vector[3]
    s_yz = 0.5 * strain_vector[4]
    s_xz = 0.5 * strain_vector[5]
    strain_tensor = numpy.array(
        [[s_xx, s_xy, s_xz], [s_xy, s_yy, s_yz], [s_yz, s_yz, s_zz]]
    )
    return strain_tensor


def write_json(filename, data_dict):
    with open(filename, "w") as fo:
        json.dump(data_dict, fo, indent=2)


def get_runtime_data(data):
    """level: "meso" or "micro"
    data: dict from runtime json
    """
    nl = data["multiscale_levels"]
    nt = data["nr_timesteps"]
    nm = data[f"{level}_nr_modes"]
    np = data[f"{level}_nr_points"] - 1  # TODO Check the +1 in points
    logger.info(f"   - detected: {nt} steps, {nm} modes, {np} points")
    return nl, nt, nm, np


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

    model = KratosMultiphysics.Model()
    parameters = KratosMultiphysics.Parameters(json.dumps(parameters_dict))
    simulation = StructuralMechanicsAnalysis(model, parameters)
    simulation.Initialize()
    modelpart = simulation._GetSolver().GetComputingModelPart()
    return model, modelpart


class Reconstruct(Common):
    def __init__(self, **kargs):
        super().__init__(**kargs)
        self.nr_voigt_comps = 6

    def element_map(self):
        """Compute auxiliar vector with the index of an element in the global vector of dofs.

        Returns:
            dict -- for each element, location of beginning in the global dof vector
            dict -- for each element, number of integration points
        """
        count = 0
        nr_voigt_comps = 6
        elem_map = {}
        nips = {}
        for element in self.modelpart.Elements:
            elem_map[element.Id] = count
            nr_ip = len(
                element.CalculateOnIntegrationPoints(
                    KratosMultiphysics.INTEGRATION_WEIGHT, self.modelpart.ProcessInfo
                )
            )
            nips[element.Id] = nr_ip
            count += nr_ip * nr_voigt_comps
        return elem_map, nips

    def get_mesh(self, rve_model):
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


    def get_material_properties(self, props):
        material_properties = {}
        material_element_list = {}
        for m in props:
            material_name = m["model_part_name"]
            logger.debug("   - loading material {}".format(material_name))
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
            for elem in self.model[material_name].Elements:
                material_element_list[material_name].append(elem.Id)
        material_elem_map = {}
        for k, v in material_element_list.items():
            for idx in v:
                material_elem_map[idx] = k
        return material_properties, material_elem_map


    def reconstruc(self, runtime_data_path):
        def initialize_runtime_micro(m, eips, nr_timesteps, nr_modes, nr_points):
            """Creates and initializes the micro runtime data file."""
            D = {}
            for x in m:
                eip = eips[x]
                e = eip[0]
                i = eip[1]
                d = {}
                d["nr_timesteps"] = nr_timesteps
                d["micro_nr_modes"] = nr_modes
                d["micro_nr_points"] = nr_points + 1
                d["micro_interpolation_parameters"] = []
                d["micro_macro_strain"] = []
                d["micro_r_value"] = []
                D[x] = d

                filename = f"micro_runtime_{e}_{i}.json"
                try:
                    os.remove(filename)
                except OSError:
                    pass

            return D

        def append_runtime_micro(D, m, eips, strain_e, c, rvalues):
            """Appends data to the micro runtime data file.
            eips is a list of tuples, containing the elemenet and the ip of the micro
            e.g. [(22, 0), (34, 7), (44, 3), ...]
            """
            for x in m:
                d = D[x]
                eip = eips[x]
                e = eip[0]
                i = eip[1]
                d["micro_interpolation_parameters"].append(c[x])  # AD_HOC FIXME
                d["micro_macro_strain"].append(list(strain_e[e][i, :]))
                d["micro_r_value"].append(rvalues[x])
                D[x] = d
            return D

        def write_runtime_micro(D, m, eips,):
            for x in m:
                e = eips[x][0]
                i = eips[x][1]
                write_json(f"micro_runtime_{e}_{i}.json", D[x])

        # Load required data meso and micro
        logger.debug("Loading runtime data {}".format(runtime_data_path))
        data = json.loads(runtime_data_path.read_text())
        levels = data["multiscale_levels"]
        nr_timesteps = data["nr_timesteps"]

        if levels == 1:
            level = "micro"
        else:
            level = "meso"
            sublevel = "micro"

        rve_interp_params = numpy.array(data[f"{level}_interpolation_parameters"])
        rve_macro_strain = numpy.array(data[f"{level}_macro_strain"])
        nr_modes = data[f"{level}_nr_modes"]
        nr_points = data[f"{level}_nr_points"] - 1
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
        r_value_correl = self.get_dataset("CORRELATION", "RVALUE", nr_modes, nr_points)
        logger.debug("Done")

        logger.debug("Loading rve data")
        rve_data = self.get_dataset("DATASET", "RVE", nr_modes, nr_points)

        logger.debug("Loading rve model")
        dset = self.get_dataset("TEMPLATE", "MODEL")
        logger.debug("Done")

        p_model = Path("model.mdpa")
        p_model.write_text(dset)
        rve_points, rve_cells = self.get_mesh(str(p_model))

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
        material_properties, material_elem_map = self.get_material_properties(
            rve_data["material_parameters"]["properties"]
        )

        ip_elem_map, nr_of_ips = self.element_map()

        # Select and initialize micros to write
        # AD-HOC. FIXME.
        if levels > 1:
            eee = [ 18, 32, 41, 4, 49, 61, 55, 68, 67, 60, 48, 10,
                    33, 61, 63, 1, 66, 4, 3, 67, 3, 4, 51, 53, 65, ]
            ppp = [ 2, 0, 0, 0, 0, 3, 3, 5, 3, 5, 7, 3, 0,
                    7, 7, 2, 4, 7, 1, 5, 3, 5, 6, 7, 2, ]
            eips = [x for x in zip(eee, ppp)]
            micro_elems = [3, 0, 2, 4, 5, ]  # elements: 4, 18, x, 41, 49, 61
            micro_nr_modes = data[f"{sublevel}_nr_modes"]
            micro_nr_points = data[f"{sublevel}_nr_points"] - 1
            mD = initialize_runtime_micro(
                micro_elems, eips, nr_timesteps, micro_nr_modes, micro_nr_points
            )

        # Open XDMF file for writing field data for each timestep
        filename = "rve_reconstructed.xdmf"
        meshio.write_points_cells(filename, rve_points, rve_cells)
        with meshio.xdmf.TimeSeriesWriter(filename) as writer:
            writer.write_points_cells(rve_points, rve_cells)
            for t in range(nr_timesteps):
                logger.info("Timestep {}".format(t))

                logger.debug("Solving fluctuant displacement")
                displacement = numpy.dot(
                    strain_correl[:, :nr_modes], rve_interp_params[t, :]
                )
                displacement = numpy.reshape(displacement, (-1, 3))

                logger.debug("Solving total displacement")
                strain_macro = rve_macro_strain[t, :]
                strain_macro_tensor = strain_voigt_to_tensor(strain_macro)
                comp = numpy.dot(strain_macro_tensor, rve_points.T)
                total_displacement = comp.T + displacement

                stress_list = self.compute_field_stress(
                    f"{level}_stress", stress_correl, data, t, nr_of_ips
                )
                strain_e, strain_h = self.compute_field_strain(
                    f"{level}_strain", strain_modes, rve_interp_params, t, nr_of_ips
                )

                ### Adding damage START
                logger.debug("Solving damage")
                damage_list = []
                r = numpy.dot(r_value_correl, data["micro_r_value"][t])
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
                    inf_yield_stress = material_properties[material_elem_map[elem_id]][
                        "inf_yield_stress"
                    ]
                    H0 = material_properties[material_elem_map[elem_id]]["H0"]
                    H1 = material_properties[material_elem_map[elem_id]]["H1"]
                    r0 = yield_stress / math.sqrt(E)
                    damage = 0
                    for r in r_in_elem[elem_id]:
                        if r < r0:
                            r = r0
                        d = 1 - q(r, E, yield_stress, inf_yield_stress, H0, H1) / r
                        damage += d / nr_ips
                    damage_list.append(damage)
                element_damage = numpy.array(damage_list).reshape(
                    (-1, 1)
                )  # formatting for meshio
                ### Adding damage END

                if levels > 1:
                    # Append micro runtime data
                    mc = data["micro_interpolation_parameters"][t]
                    mr = data["micro_r_value"][t]
                    mD = append_runtime_micro(mD, micro_elems, eips, strain_e, mc, mr)
                # Append XDMF Paraview data
                logger.debug("Writing timestep data")
                writer.write_data(
                    t,
                    point_data={
                        "DISPLACEMENT_FLUCT": numpy.reshape(displacement, (-1, 3)),
                        "DISPLACEMENT": total_displacement,
                    },
                    cell_data={
                        "STRAIN": strain_h,
                        "STRESS": stress_list,
                        "DAMAGE": element_damage,
                    },
                )
        if levels > 1:
            write_runtime_micro(mD, micro_elems, eips,)

    def compute_field_stress(self, field, stress_correl, data, t, nr_of_ips):
        logger.debug(f"Computing {field} field")
        stress = numpy.dot(stress_correl, numpy.reshape(data[field][t], (-1, 1)))
        stress_r = stress.reshape((-1, 6))
        stress_list = []
        for elem_id, nr_ips in nr_of_ips.items():
            stress_local = numpy.mean(stress_r[:nr_ips, :], axis=0)
            stress_r = stress_r[nr_ips:, :]
            stress_list.append(stress_local)
        return stress_list

    def compute_field_strain(
        self, field, strain_modes, rve_interpolation_params, t, nr_of_ips
    ):
        logger.debug(f"Computing {field} field")
        strain_global = numpy.dot(strain_modes, rve_interpolation_params[t, :])
        strain_r = strain_global.reshape((-1, 6))
        # new
        strain_e = {}
        strain_h = []
        idx = 0
        for elem_id, nr_ips in nr_of_ips.items():
            strain_e[elem_id] = strain_r[idx : idx + nr_ips, :]
            strain_h.append(numpy.mean(strain_e[elem_id], axis=0))
            idx += nr_ips
        # end new
        # strain_list = []
        # for elem_id, nr_ips in nr_of_ips.items():
        #    strain_local = numpy.mean(strain_r[:nr_ips,:], axis=0)
        #    strain_r = strain_r[nr_ips:,:]
        #    strain_list.append(strain_local)
        #    return strain_list
        return strain_e, strain_h


#####
##### wip for micro
#####

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
