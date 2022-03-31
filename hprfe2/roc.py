"""
ROC: pending description here.
"""
import logging
from pathlib import Path
import json
import numpy as np
from common import Common


logger = logging.getLogger(__name__)

###########################################################
#   Reduced Order Cubeture algorithm
###########################################################


def remove_exact_integral_energy(modes, weights):
    eps = np.finfo(float).eps  # 2.22044604925e-16
    # Total microscale volume
    total_weight = np.sum(weights)
    sqrt_total_weight = np.sqrt(total_weight)
    sqrt_weights = np.sqrt(weights)
    # Normalized exact integral
    norm_exact_integral = modes.T @ weights / total_weight
    # Matrix of modified modes (with zero integral)
    modified_modes = (modes - norm_exact_integral) * sqrt_weights.reshape(-1, 1)
    [modified_modes, bases_weights] = np.linalg.svd(
        modified_modes, full_matrices=False
    )[:2]
    # filter the reduced modified set of modes
    tolerance = np.max(modes.shape) * eps * np.max(bases_weights)
    # DEBUG
    #rank_mod_modes = sum(i > tolerance for i in bases_weights)
    rank_mod_modes = len(bases_weights)  # DEBUG
    # END DEBUG
    modified_modes = modified_modes[:, 0:rank_mod_modes]
    # Adding last row related with the sqrt of gauss integration weigths
    # and initializing the RHS vector for the optimization problem
    J = np.vstack([modified_modes.T, (sqrt_weights / sqrt_total_weight).T])
    b = np.vstack([np.zeros((modified_modes.T.shape[0], 1)), sqrt_total_weight])
    return J, b


def update_weights_inverse(H, alpha, bases_current, base_new, r):
    c = np.dot(bases_current.T, base_new)
    d = np.dot(H, c).reshape(-1, 1)
    s = np.dot(base_new.T, base_new) - np.dot(c.T, d)
    aux1 = np.hstack([H + np.outer(d, d) / s, -d / s])
    aux2 = np.hstack([np.squeeze(-d.T / s), 1 / s])
    H_new = np.vstack([aux1, aux2])
    v = np.dot(base_new.T, r) / s
    alpha = np.vstack([(alpha - d * v), v])
    return H_new, alpha


def update_inverse_hermitian(invH, neg_index):
    if neg_index == np.shape(invH)[1]:
        aux = (invH[0:-1, -1] * invH[-1, 0:-1]) / invH(-1, -1)
        invH_new = invH[:-1, :-1] - aux
    else:
        aux1 = np.hstack(
            [
                invH[:, 0:neg_index],
                invH[:, neg_index + 1 :],
                invH[:, neg_index].reshape(-1, 1),
            ]
        )
        aux2 = np.vstack(
            [aux1[0:neg_index, :], aux1[neg_index + 1 :, :], aux1[neg_index, :]]
        )
        invH_new = (
            aux2[0:-1, 0:-1] - np.outer(aux2[0:-1, -1], aux2[-1, 0:-1]) / aux2[-1, -1]
        )
    return invH_new


def multiupdate_inverse_hermitian(invH, neg_indexes):
    neg_indexes = np.sort(neg_indexes)
    for i in range(np.size(neg_indexes)):
        neg_index = neg_indexes[i] - i
        invH = update_inverse_hermitian(invH, neg_index)
    return invH


def compute_roq(Modes, weights, nGP, tol):
    J, b = remove_exact_integral_energy(Modes, weights)
    y = np.arange(len(weights))
    r = b  # residual vector, initial guess
    it = 0  # number of iterations
    mPOS = 0  # number of non-zero weights
    z = []
    Jnorm = np.sqrt(sum(np.multiply(J, J), 0))

    # point selection algorithm
    while (np.linalg.norm(r) / np.linalg.norm(b) > tol) and (mPOS <= nGP):
        # 1. Compute new point
        ObjFun = np.dot((J[:, y]).T, r)
        div = np.multiply(Jnorm[y], np.linalg.norm(r)).reshape(-1, 1)
        ObjFun = np.divide(ObjFun, div)
        s = ObjFun.argmax()
        i = y[s]
        # 2. Update alpha and H (unrestricted least squares)
        if it == 0:
            # complies with newer versions of numpy
            alpha = np.linalg.lstsq(J[:, [i]], b, rcond=None)[0]
            # alpha = np.linalg.lstsq(J[:, [i]], b)[0]
            H = 1 / np.dot((J[:, i]).T, J[:, i])
        else:
            H, alpha = update_weights_inverse(H, alpha, J[:, z], J[:, i], r)
        # 3. Move i from set y to set z
        z = (np.append(z, i)).astype(int)
        y = np.delete(y, s)
        # 4. Find possible negative weights
        if any(alpha < 0):
            logger.warning("NEGATIVE weight found")
            indexes_neg_weight = np.where(alpha <= 0.0)[0]
            y = np.append(y, (z[indexes_neg_weight]).T)
            z = np.delete(z, indexes_neg_weight)
            H = multiupdate_inverse_hermitian(H, indexes_neg_weight)
            alpha = np.dot(H, np.dot(J[:, z].T, b))

        # 6. Update the residual
        r = b - np.dot(J[:, z], alpha)
        # 7. Update mPOS and k
        mPOS = np.size(z)
        it = it + 1
        logger.debug(
            "k = {}, mPOS = {}, error = {:.2f}%".format(
                it, mPOS, np.linalg.norm(r) / np.linalg.norm(b) * 100
            )
        )
    # 6. Postprocess of points - neglecting null weights
    w = np.multiply(alpha, np.sqrt(weights[z]).reshape(-1, 1))
    logger.debug("Reduced Weights: {}".format(w.T))
    logger.debug("sum of reduced weights: {}".format(np.sum(w)))
    logger.debug("IP's index (ids starts from zero): {}".format(z))
    return w, z


###########################################################
#   End of algorithm
###########################################################


def compute_hprom_weights(ip_data, nr_roq_points, energy_modes):
    logger.info("Computing reduced set of integration points (HPROM)")
    ip_weights = ip_data[0]
    ip_lids = ip_data[1]
    elem_ids = ip_data[2]
    [w, z] = compute_roq(energy_modes, np.array(ip_weights), nr_roq_points, tol=1.0e-14)
    roq_list = []
    for x, ip_gid in enumerate(z):
        ip_lid = ip_lids[ip_gid]
        elem_id = elem_ids[ip_gid]
        roq_list.append([elem_id, ip_lid, w[x][0], ip_gid])
    return roq_list  # returns list of: element id, local IP id, IP weight, global IP id


def compute_rom_weights(ip_data):
    logger.info("Computing complete set of integration points (ROM)")
    ip_weights = ip_data[0]
    ip_lids = ip_data[1]
    elem_ids = ip_data[2]
    roq_list = []
    for ip_gid in range(len(ip_weights)):
        ip_weight = ip_weights[ip_gid]
        ip_lid = ip_lids[ip_gid]
        elem_id = elem_ids[ip_gid]
        roq_list.append([elem_id, ip_lid, ip_weight, ip_gid])
    return roq_list


def write_ip_sets(common):
    """docstring here"""

    import KratosMultiphysics
    import KratosMultiphysics.MultiscaleROMApplication
    from KratosMultiphysics.StructuralMechanicsApplication import (
        structural_mechanics_analysis,
    )

    case = common.training_path
    params = json.loads((case / "ProjectParameters_sampling.json").read_text())

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

    parameters = KratosMultiphysics.Parameters(json.dumps(params))
    model = KratosMultiphysics.Model()
    simulation = structural_mechanics_analysis.StructuralMechanicsAnalysis(
        model, parameters
    )
    simulation.Initialize()
    rve_modelpart = simulation._GetSolver().GetComputingModelPart()

    #
    # gather global model part info
    #
    ip_weights = []
    ip_lids = []
    elem_ids = []
    for elem in rve_modelpart.Elements:
        iw_list = elem.CalculateOnIntegrationPoints(
            KratosMultiphysics.INTEGRATION_WEIGHT, rve_modelpart.ProcessInfo
        )
        for ip_lid, ip_weight in enumerate(iw_list):
            ip_weights.append(ip_weight)
            ip_lids.append(ip_lid)
            elem_ids.append(elem.Id)
    ip_data = [ip_weights, ip_lids, elem_ids]
    nr_ips = len(ip_data[0])

    #
    # compute ip set
    #
    for nr_p in common.ip_subsets:
        roc_filename = common.bases_path / common.roc_fname(nr_p)
        if common.skip_calculation(roc_filename):
            logger.info("File {} exists. Skipping calculation".format(roc_filename.name))
            continue
        if "ROM" in str(nr_p):  # ROM case
            roc_list = compute_rom_weights(ip_data)
        else:  # HPROM case
            logger.info("Generating {}".format(roc_filename))
            # compute ROC list
            energy_bases = common.get_dataset("BASES", "ENERGY")[:, :nr_p]
            roc_list = compute_hprom_weights(ip_data, nr_p, energy_bases)

        with open(roc_filename, "w") as ofile:
            for list in roc_list:
                ofile.write("{} {} {} {}\n".format(list[0], list[1], list[2], list[3]))

    return


#######################################################################
# Main
#######################################################################

if __name__ == "__main__":

    import sys

    if len(sys.argv) > 1:
        C = Common(Path(sys.argv[1]))
    else:
        exit("Usage: python roc.py <root_path>")

    write_ip_sets(C)
