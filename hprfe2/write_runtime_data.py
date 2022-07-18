import os
import json
import KratosMultiphysics as KM
import KratosMultiphysics.MultiscaleROMApplication as MSR


def get_udata(model_part, element, ip):
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

    return ts, nm, np, ic, mstrain, stress, rv


class RuntimeData:
    def __init__(self, filename):
        try:
            os.remove(filename)
        except OSError:
            pass

        self.filename = filename
        self.data = {}
        self.data["nr_timesteps"] = -1
        self.data["nr_modes"] = -1
        self.data["nr_points"] = -1
        self.data["strain_coeffs"] = []
        self.data["macro_strain"] = []
        self.data["stress"] = []
        self.data["r_value"] = []

        self.data["u_nr_modes"] = -1
        self.data["u_nr_points"] = -1
        self.data["u_strain_coeffs"] = []
        self.data["u_r_value"] = []
        self.data["u_stress"] = []

        with open(filename, "w") as f:
            json.dump(self.data, f, indent=2)

    def write(self, ts, nm, np, ic, mstrain, stress, rv):

        with open(self.filename) as f:
            data = json.load(f)

        data["nr_timesteps"] = ts
        data["nr_modes"] = nm
        data["nr_points"] = np
        data["strain_coeffs"].append(ic)
        data["macro_strain"].append(mstrain)
        data["stress"].append(stress)
        data["r_value"].append(rv)

        with open(self.filename, "w") as f:
            json.dump(data, f, indent=2)

    def write_u(self, nm, np, smc, rvalue, stress):

        with open(self.filename) as f:
            data = json.load(f)

        data["u_nr_modes"] = nm
        data["u_nr_points"] = np
        data["u_strain_coeffs"].append(smc)
        data["u_r_value"].append(rvalue)
        data["u_stress"].append(stress)

        with open(self.filename, "w") as f:
            json.dump(data, f, indent=2)
