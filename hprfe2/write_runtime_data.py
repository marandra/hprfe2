import KratosMultiphysics as KM
import KratosMultiphysics.MultiscaleROMApplication as MSR
import os
import json


def write_json(filename, data_dict):
    with open(filename, "w") as fo:
        json.dump(data_dict, fo, indent=2)


def read_json(filename):
    with open(filename) as f:
        data_dict = json.load(f)
    return data_dict


def append_to_json(filename, new_data):
    data = read_json(filename)
    data["nr_timesteps"] = new_data["nr_timesteps"]
    data["nr_modes"] = new_data[f"nr_modes"]
    data["nr_points"] = new_data[f"nr_points"]
    data["interpolation_parameters"].append(new_data["interpolation_parameters"])
    data["macro_strain"].append(new_data["macro_strain"])
    data["stress"].append(new_data["stress"])
    data["r_value"].append(new_data["r_value"])
    write_json(filename, data)


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
        # self.data["multiscale_levels"] = self.multiscale_levels
        self.data["nr_timesteps"] = 0
        self.data["nr_modes"] = 0
        self.data["nr_points"] = 0
        self.data["interpolation_parameters"] = []
        self.data["macro_strain"] = []
        self.data["stress"] = []
        self.data["r_value"] = []

        write_json(filename, self.data)

    def finalize_solution_step(self):
        for elem in self.model_part.Elements:
            if elem.Id == self.element:

                ###
                self.data["nr_timesteps"] = self.model_part.ProcessInfo[KM.STEP]

                ###
                ip_data = elem.CalculateOnIntegrationPoints(
                    MSR.REDUCED_MODES_WEIGHTS_L1, self.model_part.ProcessInfo
                )
                data = [x for x in ip_data[self.ip]]
                self.data["interpolation_parameters"] = data

                ###
                self.data["nr_modes"] = len(self.data["interpolation_parameters"])

                ###
                ip_data = elem.CalculateOnIntegrationPoints(
                    KM.STRAIN, self.model_part.ProcessInfo
                )
                data = [x for x in ip_data[self.ip]]
                self.data["macro_strain"] = data

                ###
                ip_data = elem.CalculateOnIntegrationPoints(
                    KM.INTERNAL_VARIABLES, self.model_part.ProcessInfo
                )
                data = [x for x in ip_data[self.ip]]
                self.data["r_value"] = data

                ###
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
                data = [ldata[x : x + nc] for x in range(0, len(ldata), nc)]
                self.data[f"stress"] = data

                ###
                self.data["nr_points"] = len(data)

        append_to_json(self.filename, self.data)
