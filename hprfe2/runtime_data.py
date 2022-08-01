import os
import json
import numpy
import KratosMultiphysics as KM
import KratosMultiphysics.MultiscaleROMApplication as MSR


#
# Functions that process data from reconstruction
#


def write_from_reconstruction(filename, step, smc, mstrain, stress, rv):
    # if we are writing runtime data from reconstruction,
    # it is a micro already, so we don't have to worry about L2 reconstruction

    # read
    with open(filename) as f:
        data = json.load(f)

    # init
    if step == 1:
        nmodes = len(smc)
        npoints = len(stress)
        data = init_l1(data, nmodes, npoints)

    # update
    data["nr_timesteps"] = step
    data = append_l1(
        data, *get_data_from_reconstruction_l1(filename, step, smc, mstrain, stress, rv)
    )

    # write
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


def get_data_from_reconstruction_l1(filename, e, i, smc, mstrain, stress, rvalue):

    # return strain_coeffs, stress, mstrain, rv
    return smc[e][i], stress[e][i], mstrain[e][i], rvalue[e][i]


#
# Functions that process data from modelpart
#


def write_from_modelpart(filename, mp, e, i, nested=False):
    # read
    with open(filename) as f:
        data = json.load(f)

    # prepare
    step = mp.ProcessInfo[KM.STEP]

    # init
    if step == 1:
        data = init_l1(data, *prepare_from_modelpart_l1(mp, e, i))
        if nested:
            data = init_l2(data, *prepare_from_modelpart_l2(mp, e, i))

    # update
    data["nr_timesteps"] = step
    data = append_l1(data, *get_data_from_modelpart_l1(mp, e, i))
    if nested:
        data = append_l2(data, *get_data_from_modelpart_l2(mp, e, i))

    # write
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


def prepare_from_modelpart_l1(mp, e, ip):
    for elem in mp.Elements:
        if elem.Id == e:
            smc = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L1, mp.ProcessInfo
            )[ip]
            nmodes = len(smc)
            del smc

            ### volver aqui cuando optimicemos el resize de stress
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.CAUCHY_STRESS_VECTOR_L1, mp.ProcessInfo
            )[ip]
            ldata = [x for x in ip_data]
            nc = 6  # hardoced nr of comps
            # equivalent to np.resize(-1, nc) but for lists
            stress = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]
            npoints = len(stress)
            del ip_data, ldata, stress

            return nmodes, npoints


def get_data_from_modelpart_l1(model_part, element, ip):
    for elem in model_part.Elements:
        if elem.Id == element:

            ### strain modes coefficients
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L1, model_part.ProcessInfo
            )[ip]
            strain_coeffs = [x for x in ip_data]

            ### stress. L1: a vector of size npoints * ncomps
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.CAUCHY_STRESS_VECTOR_L1, model_part.ProcessInfo
            )
            ldata = [x for x in ip_data[ip]]
            nc = 6  # hardoced nr of comps
            # equivalent to np.resize(-1, nc) but for lists
            stress = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]

            ### initial strain received from macro scale
            ip_data = elem.CalculateOnIntegrationPoints(
                KM.STRAIN, model_part.ProcessInfo
            )
            macro_strain = [x for x in ip_data[ip]]

            ### r_value, i.e., internal variable of CL
            # must unpack structure:
            # [npoints, niv1, niv2, .. niv_npoints, iv0, iv1, ..., iv_n)
            ip_data = list(
                elem.CalculateOnIntegrationPoints(
                    MSR.INTERNAL_VARIABLES_L1, model_part.ProcessInfo
                )[ip]
            )
            npoints = int(ip_data.pop(0))
            niv = [int(x) for x in ip_data[:npoints]]
            del ip_data[:npoints]
            rvalue = []
            for n in niv:
                rvalue.append(ip_data[:n])
                del ip_data[:n]

            return strain_coeffs, stress, macro_strain, rvalue


def prepare_from_modelpart_l2(mp, e, ip):
    for elem in mp.Elements:
        if elem.Id == e:
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L2, mp.ProcessInfo
            )[ip]
            x = ip_data
            nr = x.Size1()
            nc = x.Size2()
            ldata = []
            for r in range(nr):
                data_i = []
                for c in range(nc):
                    data_i.append(x[r, c])
                ldata.append(data_i)
            smc = ldata
            nmodes = len(smc[0])

            ip_data = numpy.array(
                elem.CalculateOnIntegrationPoints(
                    MSR.CAUCHY_STRESS_VECTOR_L2, mp.ProcessInfo
                )[ip]
            )
            ip_data = ip_data.reshape((numpy.shape(ip_data)[0], -1, 6))
            stress = [[list(j) for j in i] for i in ip_data]
            npoints = len(stress[0])

            return nmodes, npoints


def get_data_from_modelpart_l2(model_part, element, ip):
    for elem in model_part.Elements:
        if elem.Id == element:

            ### strain modes coefficients
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L2, model_part.ProcessInfo
            )
            x = ip_data[ip]
            nr = x.Size1()
            nc = x.Size2()
            data = []
            for r in range(nr):
                data_i = []
                for c in range(nc):
                    data_i.append(x[r, c])
                data.append(data_i)
            smc = data

            ### stress
            data = numpy.array(
                elem.CalculateOnIntegrationPoints(
                    MSR.CAUCHY_STRESS_VECTOR_L2, model_part.ProcessInfo
                )[ip]
            )
            data = data.reshape((numpy.shape(data)[0], -1, 6))
            #  convert numpy 3D array to nested list
            stress = [[list(j) for j in i] for i in data]

            return smc, stress


#
# Agnostic functions
#


def append_l1(data, smc, stress, mstrain, rv):
    """read from and write to file at each timestep,
    to not loose data in case run is cancelled"""

    field = "strain_coeffs"
    if not field in data.keys():
        data[field] = []
    data[field].append(smc)

    field = "stress"
    if not field in data.keys():
        data[field] = []
    data[field].append(stress)

    field = "macro_strain"
    if not field in data.keys():
        data[field] = []
    data[field].append(mstrain)

    field = "r_value"
    if not field in data.keys():
        data[field] = []
    data[field].append(rv)

    return data


def append_l2(data, smc, stress):
    """read from and write to file at each timestep,
    to not loose data in case run is cancelled"""

    field = "u_strain_coeffs"
    if not field in data.keys():
        data[field] = []
    data[field].append(smc)

    field = "u_stress"
    if not field in data.keys():
        data[field] = []
    data[field].append(stress)

    return data


def init_l1(data, nmodes, npoints):
    data["nr_modes"] = nmodes
    data["nr_points"] = npoints
    return data


def init_l2(data, nmodes, npoints):
    data["u_nr_modes"] = nmodes
    data["u_nr_points"] = npoints
    return data


def init(filename):
    """Initialize the file and only an empty data structure"""

    try:
        os.remove(filename)
    except OSError:
        pass

    data = {}

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)
