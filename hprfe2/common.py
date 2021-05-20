# TODO:
# Check schema for input validation https://github.com/keleshev/schema
# Check fire for exposing objects to CLI https://github.com/google/python-fire

import logging
import json
from pathlib import Path
import h5py

logger = logging.getLogger(__name__)


def validate_config(default, user):
    """
    Validates and merges defaults and user configurations.
    Received configeration dictionaries
    """
    # Ideas for validation:
    # number of base modes < number of snapshots
    # number of base mode > number of requested modes

    # Check unique correlation pairs:
    # TODO

    # all keys in user must be present in default
    d_k = [dk for dk in default.keys()]
    for u_k in user.keys():
        if u_k not in d_k:
            raise SystemExit(
                'Not recognized key "{}" in user configuration. Exit.'.format(u_k)
            )

    # create context by merging default and user
    context = {**default, **user}
    return context


#######################################################################
#######################################################################


class Common:
    """
    TODO add docstrings
    """

    def __init__(self, root_path):

        config_fname = "config.json"

        try:
            config_user = json.loads((root_path / config_fname).read_text())[
                "config_data"
            ]
        except FileNotFoundError:
            # The first time there is no configuration file.
            config_user = {}

        defaults_basic = {
            # most frequently set
            "validation_dataset": [0],
            "rve_data_points": [100, 200],
            "rve_data_points_rom": False,
            "rve_data_modes": [20, 30],
            "reconstruction_pairs": [
                [20, 100],
                [30, 200],
            ],
        }
        defaults_advanced = {
            # example option
            # "rve_data_points_range_list": [[100, 250, 50], [200, 400, 100]],
            "rve_data_points_range_list": [],
            #
            "energy_name": "ENERGY_FREE",
            "energy_elastic_modes": 21,
            "energy_inelastic_modes": -1,
            "energy_svd_cutoff": 1e-4,
            "strain_name": "STRAIN_FLUCTUANT",
            "strain_elastic_modes": 6,
            "strain_inelastic_modes": -1,
            "strain_svd_cutoff": 1e-4,
            "rvalue_name": "R_VALUE",
            "rvalue_elastic_modes": 1,
            "rvalue_inelastic_modes": 30,
            "rvalue_svd_cutoff": 1e-4,
            "reuse_existing_files": True,
        }
        defaults_system = {
            # training files stuff
            "training_path": "sampling",
            "training_materials_fname": "materials.json",
            "training_model_fname": "model.mdpa",
            "case_path_pattern": "case_{}",
            "snapshots_fname": "snapshots.h5",
            # bases and roc sets files stuff
            "bases_path": "auxiliar",
            "bases_fname_pattern": "bases_{}_{}m.npy",
            "local_bases_fname_pattern": "bases_inelastic_local_{}.npy",
            "local_sv_fname_pattern": "sv_inelastic_local_{}.dat",
            "roc_fname_pattern": "roc_{}ip",
            # base + roc datasets stuff
            "datasets_path": "datasets",
            "rve_fname_pattern": "rve{}_{}m_{}ip.json",
            "correl_strain_pattern": "correlation_strain_{}m.npy",
            "correl_rvalue_pattern": "correlation_rvalue_{}m_{}ip.npy",
            # multiscale files stuff
            "multiscale_path": "validation",
        }

        # combine all levels of default options
        self.defaults_basic = defaults_basic  # keep it for initial dumping
        self.defaults_advanced = {**self.defaults_basic, **defaults_advanced}
        self.defaults = {**self.defaults_advanced, **defaults_system}
        self.config = validate_config(self.defaults, config_user)

        # file management
        self.root_path = root_path.resolve()
        self.config_file = self.root_path / config_fname
        self.training_path = self.root_path / self.config["training_path"]
        self.bases_path = self.root_path / self.config["bases_path"]
        self.datasets_path = self.root_path / self.config["datasets_path"]
        self.multiscale_path = self.root_path / self.config["multiscale_path"]
        #self.resources_path = self.root_path / self.config["resources_path"]
        self.resources_path = self.root_path / f"{self.root_path.name}.h5"

        # initialization of resources file
        self.init_dataset()

        # bases generation
        self.svd_cutoff = {}

        # self.energy_name = self.config["energy_name"]
        self.energy_elastic_modes = self.config["energy_elastic_modes"]
        self.energy_inelastic_modes = self.config["energy_inelastic_modes"]
        self.svd_cutoff[self.config["energy_name"]] = self.config["energy_svd_cutoff"]

        # self.strain_name = self.config["strain_name"]
        self.strain_elastic_modes = self.config["strain_elastic_modes"]
        self.strain_inelastic_modes = self.config["strain_inelastic_modes"]
        self.svd_cutoff[self.config["strain_name"]] = self.config["strain_svd_cutoff"]

        # self.rvalue_name = self.config["rvalue_name"]
        self.rvalue_elastic_modes = self.config["rvalue_elastic_modes"]
        self.rvalue_inelastic_modes = self.config["rvalue_inelastic_modes"]
        self.svd_cutoff[self.config["rvalue_name"]] = self.config["rvalue_svd_cutoff"]

        # points
        self.ip_subsets = [x for x in self.config["rve_data_points"]]
        for r in self.config["rve_data_points_range_list"]:  # unpack list of "ranges"
            for i in range(*r):
                self.ip_subsets.append(i)
        self.ip_subsets = list(set(self.ip_subsets))
        self.ip_subsets.sort()
        if self.config["rve_data_points_rom"]:
            self.ip_subsets.append("ROM")

        self.roc_fname_pattern = self.config["roc_fname_pattern"]

        # modes
        # self.reduced_nr_modes = self.config["rve_data_modes"]

        self.materials_fname = self.training_path / Path(
            self.config["training_materials_fname"]
        )
        self.rve_fname_pattern = self.config["rve_fname_pattern"]

    def dump_config(self, fname, mode=""):
        """
        Writes configuration file.
        Mode:
        'defaults_basic': defaults for common parameters (used for initial config files)
        'defaults': defaults values
        any other value dumps current config values.
        """
        if mode == "defaults_basic":
            Path(fname).write_text(
                json.dumps({"config_data": self.defaults_basic}, indent=2)
            )
        elif mode == "defaults_advanced":
            Path(fname).write_text(
                json.dumps({"config_data": self.defaults_advanced}, indent=2)
            )
        else:
            Path(fname).write_text(json.dumps({"config_data": self.config}, indent=2))

    def case_name(self, c_id):
        """
        Returns case name with corresponging leading zeros
        e.g. if nr_cases:100, id:00..99, nr_id: 2 -> case_01..case_99
        """
        # strain_set = self.parse_training_strain_set()
        # nr_cases = len(strain_set)
        nr_cases = 999  # TODO: hardcoded
        len_id = len(str(nr_cases - 1))  # size of the case number string
        case_id = "{:0{}d}".format(c_id, len_id)
        case_name = self.config["case_path_pattern"].format(case_id)
        return case_name

    def roc_fname(self, points):
        """
        docstrings here
        """
        return self.roc_fname_pattern.format(points)

    def rve_fname(self, order, modes, points):
        """
        docstrings here
        """
        return self.rve_fname_pattern.format(order, modes, points)

    def skip_calculation(self, fname):
        """
        Generates a list of files following filename pattern.
        Length of list is used as flag (False if empty, True otherwise)
        """
        fpath = Path.cwd() / fname  # converts filename to absolute Path
        flag_exists = len([f for f in fpath.parent.glob(fpath.name)])
        flag_reuse = self.config["reuse_existing_files"]
        return flag_exists and flag_reuse

    def get_bases_fname(self, field):
        """
        docstrings here
        """
        filename = self.config["bases_fname_pattern"].format(field, "*")
        fpath = self.bases_path / filename
        files = [f for f in fpath.parent.glob(fpath.name)]
        if len(files) == 0:
            return None
        if len(files) > 1:
            logger.warning(
                "More than one {} bases file detected. "
                "Picking first in the list: {}".format(field, files[0].name)
            )
        return files[0]


    def init_dataset(self):
        # Check valid keys
        #valid_keys = {}
        #valid_keys["BASES"] = ["STRAIN", "ENERGY", "RVALUE", "TEST"]
        #valid_keys["CORRELATION"] = ["STRAIN", "RVALUE"]
        #valid_keys["DATASET"] = ["RVE"]
        #valid_keys["TEMPLATE"] = ["MODEL", "MATERIALS", "PARAMETERS", "MAIN", "RESOURCES"]
        #if group not in valid_keys.keys():
        #   logger.error(f"Invalid group name {group}")
        #   exit()
        #if dataset not in valid_keys[group]:
        #   logger.error(f"Invalid dataset name {dataset}")
        #   exit()
        with h5py.File(self.resources_path, "a") as f:
            for group in ["BASES", "CORRELATION", "DATASET", "TEMPLATE"]:
                if group in f.keys():
                    continue
                f.create_group(group)

    def name_dataset(self, dataset, nmodes=None, npoints=None):
        # Set name
        name = f"{dataset}"
        if nmodes:
            name += f"_{nmodes}m"
        if nmodes and npoints:
            name += f"_{npoints}p"
        return name

    def get_dataset(self, group, dataset, nmodes=None, npoints=None):
        """Get dataset from database h5 file"""
        dsname = self.name_dataset(dataset, nmodes, npoints)
        with h5py.File(self.resources_path, "r") as f:
            if group in ["BASES", "CORRELATION"]:  # numpy
                return f[group][dsname][()]
            elif group in ["DATASET"]:  # json
                return json.loads(f[group][dsname].asstr()[()])
            else:  # text
                return f[group][dsname].asstr()[()]

    def set_dataset(self, data, group, dataset, nmodes=None, npoints=None):
        """Create dataset in database h5 file. Only valid groups and datasets"""
        dsname = self.name_dataset(dataset, nmodes, npoints)
        with h5py.File(self.resources_path, "a") as f:
            f.create_dataset(f"{group}/{dsname}", data=data)

    def has_dataset(self, group, dataset, nmodes=None, npoints=None):
        dsname = self.name_dataset(dataset, nmodes, npoints)
        with h5py.File(self.resources_path, "a") as f:
            return dsname in f[group]

#####################################################################
# main
#####################################################################

# TODO: move this to tests

if __name__ == "__main__":
    C = Common(root_path=Path("."))
    print("Test:")
    print(C.roc_fname("1"))
    print(C.roc_fname("100"))
    print(C.roc_fname("1000"))
    print(C.roc_fname(1000))
    print(C.roc_fname("ROM"))
    print(C.rve_fname(1, 20, "1"))
    print(C.rve_fname(1, 20, "100"))
    print(C.rve_fname(1, 200, "1000"))
    print(C.rve_fname(1, 200, 1000))
    print(C.rve_fname(1, "2000", "ROM"))
    print(C.ip_subsets)
