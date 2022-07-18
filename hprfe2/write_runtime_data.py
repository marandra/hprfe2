import KratosMultiphysics as KM
import KratosMultiphysics.MultiscaleROMApplication as MSR
import os
import json


def append_to_json(filename, ts, nm, np, ic, mstrain, stress, rv):
    with open(filename) as f:
        data = json.load(f)

    data["nr_timesteps"] = ts
    data["nr_modes"] = nm
    data["nr_points"] = np
    data["interpolation_parameters"].append(ic)
    data["macro_strain"].append(mstrain)
    data["stress"].append(stress)
    data["r_value"].append(rv)

    with open(filename, "w") as f:
        json.dump(data, f, indent=2)


class RuntimeData:
    def __init__(self, filename, model_part, elem, ip):
        try:
            os.remove(filename)
        except OSError:
            pass

        self.filename = filename
        self.model_part = model_part
        self.element = elem
        self.ip = ip

        self.data = {}
        self.data["nr_timesteps"] = 0
        self.data["nr_modes"] = 0
        self.data["nr_points"] = 0
        self.data["interpolation_parameters"] = []
        self.data["macro_strain"] = []
        self.data["stress"] = []
        self.data["r_value"] = []

        with open(filename, "w") as f:
            json.dump(self.data, f, indent=2)

    def finalize_solution_step(self):
        for elem in self.model_part.Elements:
            if elem.Id == self.element:

                ### timestep
                ts = self.model_part.ProcessInfo[KM.STEP]

                ### interpolation coef of strain modes
                ip_data = elem.CalculateOnIntegrationPoints(
                    MSR.REDUCED_MODES_WEIGHTS_L1, self.model_part.ProcessInfo
                )
                ic = [x for x in ip_data[self.ip]]

                ### number of modes
                nm = len(ic)

                ### initial strain received from macro scale
                ip_data = elem.CalculateOnIntegrationPoints(
                    KM.STRAIN, self.model_part.ProcessInfo
                )
                mstrain = [x for x in ip_data[self.ip]]

                ### r_value, i.e., internal variable of CL
                ip_data = elem.CalculateOnIntegrationPoints(
                    KM.INTERNAL_VARIABLES, self.model_part.ProcessInfo
                )
                rv = [x for x in ip_data[self.ip]]

                ### stress
                ip_data = elem.CalculateOnIntegrationPoints(
                    MSR.CAUCHY_STRESS_VECTOR_L1, self.model_part.ProcessInfo
                )
                # x = ip_data[self.ip]
                # nr = x.Size1()
                # nc = x.Size2()
                # data = []
                # for r in range(nr):
                #    data_i = []
                #    for c in range(nc):
                #        data_i.append(x[r, c])
                #    data.append(data_i)
                ldata = [x for x in ip_data[self.ip]]
                nc = 6  # hardoced nr of comps
                stress = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]

                ### number of points
                np = len(stress)

        append_to_json(self.filename, ts, nm, np, ic, mstrain, stress, rv)
