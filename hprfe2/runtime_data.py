import os
import json
import KratosMultiphysics as KM
import KratosMultiphysics.MultiscaleROMApplication as MSR


def get_udata(model_part, element, ip):
def get_data_l2(model_part, element, ip):
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
            ### number of points
            np = len(smc)
            ### number of modes
            nm = len(smc[0])

            ## r_value, i.e., internal variable of CL
            ## ther eno L2, because RL returns all iv in one vector
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.INTERNAL_VARIABLES_L1, model_part.ProcessInfo
            )
            data = [x for x in ip_data[ip]]
            rvalue = data

            ### stress
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.CAUCHY_STRESS_VECTOR_L2, model_part.ProcessInfo
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
            #ldata = [x for x in ip_data[ip]]
            #nc = 6  # hardoced nr of comps
            #stress = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]
            stress = data

    return nm, np, smc, rvalue, stress


def get_data(model_part, element, ip):
def get_data_l1(model_part, element, ip):
    for elem in model_part.Elements:
        if elem.Id == element:

            ### timestep
            ts = model_part.ProcessInfo[KM.STEP]

            ### strain modes coefficients
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.REDUCED_MODES_WEIGHTS_L1, model_part.ProcessInfo
            )
            ic = [x for x in ip_data[ip]]

            ### number of modes
            nm = len(ic)

            ### initial strain received from macro scale
            ip_data = elem.CalculateOnIntegrationPoints(
                KM.STRAIN, model_part.ProcessInfo
            )
            mstrain = [x for x in ip_data[ip]]

            ### r_value, i.e., internal variable of CL
            ip_data = elem.CalculateOnIntegrationPoints(
                KM.INTERNAL_VARIABLES, model_part.ProcessInfo
            )
            rv = [x for x in ip_data[ip]]

            ### stress
            ip_data = elem.CalculateOnIntegrationPoints(
                MSR.CAUCHY_STRESS_VECTOR_L1, model_part.ProcessInfo
            )
            # x = ip_data[ip]
            # nr = x.Size1()
            # nc = x.Size2()
            # data = []
            # for r in range(nr):
            #    data_i = []
            #    for c in range(nc):
            #        data_i.append(x[r, c])
            #    data.append(data_i)
            ldata = [x for x in ip_data[ip]]
            nc = 6  # hardoced nr of comps
            stress = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]

            ### number of points
            np = len(stress)

def append_l1(data, nm, np, smc, mstrain, stress, rv):
    """read from and write to file at each timestep,
    to not loose data in case run is cancelled"""

    data["nr_modes"] = nm
    data["nr_points"] = np
    data["strain_coeffs"].append(smc)
    data["macro_strain"].append(mstrain)
    data["stress"].append(stress)
    data["r_value"].append(rv)
    return data


def append_l2(data, nm, np, smc, rvalue, stress):
    """read from and write to file at each timestep,
    to not loose data in case run is cancelled"""

    data["u_nr_modes"] = len(smc[0])
    data["u_nr_points"] = len(smc[0])
    data["u_strain_coeffs"].append(smc)
    data["u_r_value"].append(rvalue)
    data["u_stress"].append(stress)
    return data



def write(filename, modelpart, e, i, nested):

    # read
    with open(filename) as f:
        data = json.load(f)

    step = modelpart.ProcessInfo[KM.STEP]

    # init
    if step == 1:
        data = init_l1(data, modelpart, e, i)
        if nested:
            data = init_l2(data, modelpart, e, i)

    # update
    data["nr_timesteps"] = step
    data = append_l1(data, *get_data_l1(modelpart, e, i))
    if nested:
        data = append_l2(data, *get_data_l2(modelpart, e, i))

    # write
    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


def init_l1(data, mp, e, ip):

    for elem in mp.Elements:
        if elem.Id == e:
            smc = elem.CalculateOnIntegrationPoints(
                    MSR.REDUCED_MODES_WEIGHTS_L1,
                    mp.ProcessInfo
                    )[ip]
            nmodes = len(smc)
        break

    data["nr_modes"] = nmodes
    data["nr_points"] = -1

    #data["strain_coeffs"] = []
    #data["macro_strain"] = []
    #data["stress"] = []
    #data["r_value"] = []
    
    return data


def init_l2(data, modelpart, e, i):

    data["u_nr_modes"] = -1
    data["u_nr_points"] = -1

    #data["u_strain_coeffs"] = []
    #data["u_r_value"] = []
    #data["u_stress"] = []
    
    return data


def init(filename, modelpart, e, i, nested):

    try:
        os.remove(filename)
    except OSError:
        pass

    data = {}

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)
