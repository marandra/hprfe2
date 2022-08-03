import os
import json
import numpy
import KratosMultiphysics as KM
import KratosMultiphysics.MultiscaleROMApplication as MSR

# Standarize nomeclature

NSTEPS = "nr_timesteps"  # step
NMODES = "nr_modes"  # nmodes
NPOINTS = "nr_points"  # npoints
CSTRAIN = "strain_coeffs"  # cstrain
STRESS = "stress"  # stress
RVALUE = "rvalue"  # rvalue
MSTRAIN = "macro_strain"  # mstrain
UNMODES = "u_" + NMODES
UNPOINTS = "u_" + NPOINTS
UCSTRAIN = "u_" + CSTRAIN
USTRESS = "u_" + STRESS
URVALUE = "u_" + RVALUE
UMSTRAIN = "u_" + MSTRAIN

#
# Getters
#


def get_nsteps(data):
    return data[NSTEPS]


def get_nmodes(data):
    return data[NMODES]


def get_npoints(data):
    return data[NPOINTS] - 1


def get_cstrain(data, step=None):
    if step:
        return data[CSTRAIN][step - 1]
    return data[CSTRAIN]


def get_rvalue(data, step=None):
    if step:
        return data[RVALUE][step - 1]
    return data[RVALUE]


def get_mstrain(data, step=None):
    if step:
        return data[MSTRAIN][step - 1]
    return data[MSTRAIN]

def udata(data):
    return UNPOINTS in data.keys()
#
# Functions that process data from reconstruction
#


def write_from_reconstruction(fname, rtdata, mstrain, step, idx):
    # - if we are writing runtime data from reconstruction,
    # it is a micro already, so we don't have to worry about L2 reconstruction
    # - step starts from 1

    # read
    with open(fname) as f:
        data = json.load(f)

    # init
    if step == 1:
        data = init_l1(data, *prepare_from_reconstruction_l1(rtdata, idx))

    # update
    data[NSTEPS] = step
    data = append_l1(data, *get_data_from_reconstruction_l1(rtdata, step, idx), mstrain)

    # write
    with open(fname, "w") as f:
        json.dump(data, f, indent=2)


def prepare_from_reconstruction_l1(rtdata, idx):
    nmodes = len(rtdata[UCSTRAIN][0][idx])
    npoints = len(rtdata[USTRESS][0][idx])
    return nmodes, npoints


def get_data_from_reconstruction_l1(rtdata, step, idx):
    cstrain = rtdata[UCSTRAIN][step - 1][idx]
    stress = rtdata[USTRESS][step - 1][idx]
    rvalue = rtdata[RVALUE][step - 1][idx]
    return cstrain, stress, rvalue


#
# Functions that process data from modelpart
#


def write_from_modelpart(fname, mp, e, ip, nested=False):
    # read
    with open(fname) as f:
        data = json.load(f)

    # prepare
    step = mp.ProcessInfo[KM.STEP]

    # init
    if step == 1:
        data = init_l1(data, *prepare_from_modelpart_l1(mp, e, ip))
        if nested:
            data = init_l2(data, *prepare_from_modelpart_l2(mp, e, ip))

    # update
    data[NSTEPS] = step
    data = append_l1(data, *get_data_from_modelpart_l1(mp, e, ip))
    if nested:
        data = append_l2(data, *get_data_from_modelpart_l2(mp, e, ip))

    # write
    with open(fname, "w") as f:
        json.dump(data, f, indent=2)


def prepare_from_modelpart_l1(mp, e, ip):
    for elem in mp.Elements:
        if elem.Id == e:
            cstrain = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L1, mp.ProcessInfo
            )[ip]
            nmodes = len(cstrain)
            del cstrain

            ### TODO: volver aqui cuando optimicemos el resize de stress
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


def get_data_from_modelpart_l1(mp, e, ip):
    for elem in mp.Elements:
        if elem.Id == e:

            ### strain modes coefficients
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L1, mp.ProcessInfo
            )[ip]
            cstrain = [x for x in ip_data]

            ### stress. L1: a vector of size npoints * ncomps
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.CAUCHY_STRESS_VECTOR_L1, mp.ProcessInfo
            )
            ldata = [x for x in ip_data[ip]]
            nc = 6  # hardoced nr of comps
            # equivalent to np.resize(-1, nc) but for lists
            stress = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]

            ### initial strain received from macro scale
            ip_data = elem.CalculateOnIntegrationPoints(KM.STRAIN, mp.ProcessInfo)
            mstrain = [x for x in ip_data[ip]]

            ### r_value, i.e., internal variable of CL
            # must unpack structure:
            # [npoints, niv1, niv2, .. niv_npoints, iv0, iv1, ..., iv_n)
            ip_data = list(
                elem.CalculateOnIntegrationPoints(
                    MSR.INTERNAL_VARIABLES_L1, mp.ProcessInfo
                )[ip]
            )
            npoints = int(ip_data.pop(0))
            niv = [int(x) for x in ip_data[:npoints]]
            del ip_data[:npoints]
            rvalue = []
            for n in niv:
                rvalue.append(ip_data[:n])
                del ip_data[:n]

            return cstrain, stress, rvalue, mstrain


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
            cstrain = ldata
            nmodes = len(cstrain[0])

            ip_data = numpy.array(
                elem.CalculateOnIntegrationPoints(
                    MSR.CAUCHY_STRESS_VECTOR_L2, mp.ProcessInfo
                )[ip]
            )
            ip_data = ip_data.reshape((numpy.shape(ip_data)[0], -1, 6))
            stress = [[list(j) for j in i] for i in ip_data]
            npoints = len(stress[0])

            return nmodes, npoints


def get_data_from_modelpart_l2(mp, e, ip):
    for elem in mp.Elements:
        if elem.Id == e:

            ### strain modes coefficients
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L2, mp.ProcessInfo
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
            cstrain = data

            ### stress
            data = numpy.array(
                elem.CalculateOnIntegrationPoints(
                    MSR.CAUCHY_STRESS_VECTOR_L2, mp.ProcessInfo
                )[ip]
            )
            data = data.reshape((numpy.shape(data)[0], -1, 6))
            #  convert numpy 3D array to nested list
            stress = [[list(j) for j in i] for i in data]

            return cstrain, stress


#
# Agnostic functions
#


def append_l1(data, cstrain, stress, rvalue, mstrain):
    """read from and write to file at each timestep,
    to not loose data in case run is cancelled"""

    field = CSTRAIN
    if not field in data.keys():
        data[field] = []
    data[field].append(cstrain)

    field = STRESS
    if not field in data.keys():
        data[field] = []
    data[field].append(stress)

    field = RVALUE
    if not field in data.keys():
        data[field] = []
    data[field].append(rvalue)

    field = MSTRAIN
    if not field in data.keys():
        data[field] = []
    data[field].append(mstrain)

    return data


def append_l2(data, cstrain, stress):
    """read from and write to file at each timestep,
    to not loose data in case run is cancelled"""

    field = UCSTRAIN
    if not field in data.keys():
        data[field] = []
    data[field].append(cstrain)

    field = USTRESS
    if not field in data.keys():
        data[field] = []
    data[field].append(stress)

    return data


def init_l1(data, nmodes, npoints):
    data[NMODES] = nmodes
    data[NPOINTS] = npoints
    return data


def init_l2(data, nmodes, npoints):
    data[UNMODES] = nmodes
    data[UNPOINTS] = npoints
    return data


def init(fname):
    """Initialize the file and only an empty data structure"""

    try:
        os.remove(fname)
    except OSError:
        pass

    data = {}

    with open(fname, "w") as f:
        json.dump(data, f, indent=2)
