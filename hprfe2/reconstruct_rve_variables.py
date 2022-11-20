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
from pathlib import Path
import logging
import math
import json
import numpy as np

from docopt import docopt
import meshio

import KratosMultiphysics as KM
from KratosMultiphysics.StructuralMechanicsApplication.structural_mechanics_analysis import (
    StructuralMechanicsAnalysis,
)

from common import Common
import runtime_data as rtd

np.set_printoptions(
    linewidth=120,
    suppress=True,
)


#
# Functions for DAMAGE reconstruction
#
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
    material = {}
    elem_list = {}
    for m in props:
        name = m["model_part_name"]
        logger.debug("   - loading material {}".format(name))
        logger.debug(m["Material"]["Variables"])
        material[name] = {}
        ym = m["Material"]["Variables"]["YOUNG_MODULUS"]
        nu = m["Material"]["Variables"]["POISSON_RATIO"]
        y = m["Material"]["Variables"]["STRESS_LIMITS"][0]
        iy = m["Material"]["Variables"]["STRESS_LIMITS"][1]
        h0 = m["Material"]["Variables"]["HARDENING_PARAMETERS"][0]
        h1 = m["Material"]["Variables"]["HARDENING_PARAMETERS"][1]
        material[name]["E"] = ym
        material[name]["nu"] = nu
        material[name]["yield_stress"] = y
        material[name]["inf_yield_stress"] = iy
        material[name]["H0"] = h0
        material[name]["H1"] = h1
        material[name]["C"] = compute_elastic_tensor(ym, nu)

        elem_list[name] = []
        for elem in model[name].Elements:
            elem_list[name].append(elem.Id)
    elem_map = {}
    for k, v in elem_list.items():
        for idx in v:
            elem_map[idx] = k
    return material, elem_map


def compute_damage(rtd, data, t, rvalue_correl, ips_per_elem, material, elem_map):
    def q(r, e, y, iy, h0, h1):
        r0 = y / math.sqrt(e)
        q0 = r0  # strain_variable_init
        q1 = iy / math.sqrt(e)  # stress_variable_inf
        r1 = r0 + (q1 - q0) / h0
        if r < r0:
            return q0
        if r >= r0 and r < r1:
            return q0 + h0 * (r - r0)
        # Case r >= r1:
        return q1 + h1 * (r - r1)

    rvalue = [x[0] for x in rtd.get_rvalue(data, t + 1)]  # shape (n, 1) -> (n,)
    rvalue_global = np.dot(rvalue_correl, rvalue)
    rvalue_in_ips = {}
    for e, nip in ips_per_elem.items():
        rvalue_in_ips[e] = rvalue_global[:nip]
        rvalue_global = rvalue_global[nip:]
    damage = []
    for e, nip in ips_per_elem.items():
        # C = material[elem_map[e]]["C"]
        # nu = material[elem_map[e]]["nu"]
        ym = material[elem_map[e]]["E"]
        y = material[elem_map[e]]["yield_stress"]
        iy = material[elem_map[e]]["inf_yield_stress"]
        h0 = material[elem_map[e]]["H0"]
        h1 = material[elem_map[e]]["H1"]
        r0 = y / math.sqrt(ym)
        d_elem = 0
        for r in rvalue_in_ips[e]:
            if r < r0:
                r = r0
            d_ip = 1 - q(r, ym, y, iy, h0, h1) / r  # damage at ip
            d_elem += d_ip / nip  # damage homogenized in elem
        damage.append(d_elem)
    return np.array(damage).reshape((-1, 1))  # formatting for meshio

def compute_stress_from_damage(damage_h, strain_h, ips_per_elem, material, elem_map):
    stress_elem = []
    for e, _ in ips_per_elem.items():
        C = material[elem_map[e]]["C"]
        s = strain_h[e-1].reshape((6, 1))
        stress = (1 - damage_h[e-6]) * np.dot(C, s)
        stress_elem.append(stress)
    return np.array(stress_elem).reshape((-1, 6))  # formatting for meshio


#
# End functions for DAMAGE reconstruction
#


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
        element_type = cell_block.type
        # if "hexa" in element_type or "wedge" in element_type:
        if "line8" in element_type:
            rve_cells.append(meshio.CellBlock("hexahedron", cell_block.data))
        if "line6" in element_type:
            rve_cells.append(meshio.CellBlock("wedge", cell_block.data))
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


def ei_to_reconstr(rtdata_path):

    # Begin ADHOC. FIXME
    # Select and initialize micros to write
    uei_pairs = [
        (0, 18, 2),
        (1, 32, 0),
        (2, 41, 0),
        (3, 4, 0),
        (4, 49, 0),
        (5, 61, 3),
        (6, 55, 3),
        (7, 68, 5),
        (8, 67, 3),
        (9, 60, 5),
        (10, 48, 7),
        (11, 10, 3),
        (12, 33, 0),
        (13, 61, 7),
        (14, 63, 7),
        (15, 1, 2),
        (16, 66, 4),
        (17, 4, 7),
        (18, 3, 1),
        (19, 67, 5),
        (20, 3, 3),
        (21, 4, 5),
        (22, 51, 6),
        (23, 53, 7),
        (24, 65, 2),
    ]
    uelems = [3, 0, 2, 4, 5]  # elements: 4, 18, x, 41, 49, 61
    aux = [uei_pairs[i] for i in uelems]
    # for 35p
    aux = [(0, 1, 7), (1, 13, 7), (2, 28, 0), (3, 38, 6), (4, 53, 6), (5, 61, 5)]
    # End ADHOC

    points = []
    for p in aux:
        points.append(
            {
                "fname": f"{rtdata_path.stem}-u_{p[1]}_{p[2]}.json",
                "idx": p[0],
                "ee": p[1],
                "ii": p[2],
            }
        )
    return points


class Reconstruct(Common):
    def __init__(self, rtdata_path, **kargs):
        super().__init__(**kargs)
        self.rtdata_path = rtdata_path
        self.nr_voigt_comps = 6
        self.reconstruct_micro = False
        self.filename = f"{rtdata_path.stem}_reconstr.xdmf"

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

    def gather_rtdata(self, rtdata_path):
        data = json.loads(rtdata_path.read_text())
        nsteps = rtd.get_nsteps(data)
        nmodes = rtd.get_nmodes(data)
        npoints = rtd.get_npoints(data)
        cstrain = np.array(rtd.get_cstrain(data))
        stress = rtd.get_mstrain(data)
        rvalue = rtd.get_rvalue(data)
        mstrain = np.array(rtd.get_mstrain(data))
        # TODO: Here we check if it is a meso or a micro, move it to other place
        self.reconstruct_micro = True if rtd.udata(data) else False
        return nsteps, nmodes, npoints, cstrain, stress, rvalue, mstrain

    def gather_datasets(self, nmodes, npoints):
        strain_modes = self.get_dataset("BASES", "STRAIN")[:, :nmodes]
        strain_correl = self.get_dataset("CORRELATION", "STRAIN", nmodes)
        stress_correl = self.get_dataset("CORRELATION", "STRESS", nmodes, npoints)
        self.skip_damage_reconstruction = False
        r_value_correl = None
        try:
            r_value_correl = self.get_dataset("CORRELATION", "RVALUE", nmodes, npoints)
        except KeyError:
            logger.warning(
                "RVALUE correlation matrix not present. Skipping DAMAGE reconstruction"
            )
            self.skip_damage_reconstruction = True

        rve_data = self.get_dataset("DATASET", "RVE", nmodes, npoints)
        model = self.get_dataset("TEMPLATE", "MODEL")
        return (
            strain_modes,
            strain_correl,
            stress_correl,
            r_value_correl,
            rve_data,
            model,
        )

    def gather_model(self, model):
        p_model = Path("model.mdpa")
        p_model.write_text(model)
        rve_points, rve_cells = get_mesh(str(p_model))

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
        return rve_points, rve_cells

    def compute_field_stress(self, field, stress_correl, data, t, nr_of_ips):
        logger.debug(f" - Computing {field} field")
        # Este es el stress en cada punto de gauss meso
        stress_global = np.dot(stress_correl, np.reshape(data[field][t], (-1, 1)))
        stress_r = stress_global.reshape((-1, 6))
        stress_e = {}
        # Este es el stress homogenizado en cada elemento meso
        stress_h = []
        idx = 0
        for elem_id, nr_ips in nr_of_ips.items():
            stress_e[elem_id] = stress_r[idx : idx + nr_ips, :]
            stress_h.append(np.mean(stress_e[elem_id], axis=0))
            idx += nr_ips
        return stress_e, np.array(stress_h).reshape(-1, 6)

    def compute_field_strain(
        self, field, strain_modes, rve_interpolation_params, t, nr_of_ips
    ):
        logger.debug(f" - Computing {field} field")
        strain_global = np.dot(strain_modes, rve_interpolation_params[t, :])
        strain_r = strain_global.reshape((-1, 6))
        strain_e = {}
        strain_h = []
        idx = 0
        for elem_id, nr_ips in nr_of_ips.items():
            strain_e[elem_id] = strain_r[idx : idx + nr_ips, :]
            strain_h.append(np.mean(strain_e[elem_id], axis=0))
            idx += nr_ips
        return strain_e, np.array(strain_h).reshape((-1, 6))

    def reconstruct(self):
        rtdata_path = self.rtdata_path
        logger.debug(f"Loading runtime data {rtdata_path}")
        nsteps, nmodes, npoints, cstrain, stress, rvalue, mstrain = self.gather_rtdata(
            rtdata_path
        )

        logger.debug("Loading databases")
        (
            strain_modes,
            strain_correl,
            stress_correl,
            r_value_correl,
            rve_data,
            model,
        ) = self.gather_datasets(nmodes, npoints)

        logger.debug("Loading rve model and materials")
        rve_points, rve_cells = self.gather_model(model)

        if not self.skip_damage_reconstruction:
            material_properties, material_elem_map = get_material_properties(
                self.model, rve_data["material_parameters"]["properties"]
            )

        ip_elem_map, nr_of_ips = self.element_map()

        # Generate micro runtime data if required data present
        if self.reconstruct_micro:
            for d in ei_to_reconstr(self.rtdata_path):
                rtd.init(d["fname"])

        # Open XDMF file for writing field data for each timestep
        meshio.write_points_cells(self.filename, rve_points, rve_cells)
        with meshio.xdmf.TimeSeriesWriter(self.filename) as writer:
            writer.write_points_cells(rve_points, rve_cells)
            for t in range(nsteps):
                logger.info("Timestep {}".format(t))

                logger.debug(" - Solving fluctuant displacement")
                displacement = np.dot(strain_correl[:, :nmodes], cstrain[t, :])
                displacement = np.reshape(displacement, (-1, 3))

                logger.debug(" - Solving total displacement")
                strain_macro = mstrain[t, :]
                strain_macro_tensor = strain_voigt_to_tensor(strain_macro)
                comp = np.dot(strain_macro_tensor, rve_points.T)
                total_displacement = comp.T + displacement

                data = json.loads(rtdata_path.read_text())
                stress_e, stress_h = self.compute_field_stress(
                    f"stress", stress_correl, data, t, nr_of_ips
                )
                strain_e, strain_h = self.compute_field_strain(
                    f"strain", strain_modes, cstrain, t, nr_of_ips
                )

                if not self.skip_damage_reconstruction:
                    logger.debug(" - Solving damage")
                    data = json.loads(rtdata_path.read_text())
                    damage_h = compute_damage(
                        rtd,
                        data,
                        t,
                        r_value_correl,
                        nr_of_ips,
                        material_properties,
                        material_elem_map,
                    )
                    stressd_h = compute_stress_from_damage(damage_h, strain_h,
                        nr_of_ips,
                        material_properties,
                        material_elem_map,
                )

                if self.reconstruct_micro:
                    logger.debug(" - Writing micro runtime data")
                    data = json.loads(rtdata_path.read_text())
                    for d in ei_to_reconstr(self.rtdata_path):
                        fname = d["fname"]
                        idx = d["idx"]
                        ee = d["ee"]
                        ii = d["ii"]
                        rtd.write_from_reconstruction(
                            fname, data, list(strain_e[ee][ii]), t + 1, idx
                        )

                # Append XDMF Paraview data
                logger.debug(" - Writing timestep data")
                point_data = {}
                cell_data = {}
                point_data["DISPLACEMENT_FLUCT"] = np.reshape(displacement, (-1, 3))
                point_data["DISPLACEMENT"] = total_displacement
                workaround_flag = (
                    False  # FIXME. Workaround for paraview 5.10, ok for <=5.9
                )
                if not workaround_flag:
                    cell_data["STRAIN"] = strain_h
                    cell_data["STRESS"] = stress_h
                else:
                    cell_data["STRAIN_MAGNITUDE"] = np.linalg.norm(
                        strain_h, axis=1
                    ).reshape((-1, 1))
                    cell_data["STRAIN_XX"] = strain_h[:, 0].reshape((-1, 1))
                    cell_data["STRAIN_YY"] = strain_h[:, 1].reshape((-1, 1))
                    cell_data["STRAIN_ZZ"] = strain_h[:, 2].reshape((-1, 1))
                    cell_data["STRAIN_XY"] = strain_h[:, 3].reshape((-1, 1))
                    cell_data["STRAIN_YZ"] = strain_h[:, 4].reshape((-1, 1))
                    cell_data["STRAIN_XZ"] = strain_h[:, 5].reshape((-1, 1))
                    cell_data["STRESS_MAGNITUDE"] = np.linalg.norm(
                        stress_h, axis=1
                    ).reshape((-1, 1))
                    cell_data["STRESS_XX"] = stress_h[:, 0].reshape((-1, 1))
                    cell_data["STRESS_YY"] = stress_h[:, 1].reshape((-1, 1))
                    cell_data["STRESS_ZZ"] = stress_h[:, 2].reshape((-1, 1))
                    cell_data["STRESS_XY"] = stress_h[:, 3].reshape((-1, 1))
                    cell_data["STRESS_YZ"] = stress_h[:, 4].reshape((-1, 1))
                    cell_data["STRESS_XZ"] = stress_h[:, 5].reshape((-1, 1))
                if not self.skip_damage_reconstruction:
                    cell_data["DAMAGE"] = damage_h
                    cell_data["STRESS_D"] = stressd_h
                writer.write_data(t, point_data=point_data, cell_data=cell_data)


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

    RECONST = Reconstruct(Path(ARGS["<runtime_data>"]), root_path=Path(ARGS["<root>"]))
    RECONST.reconstruct()
