# TODO:
# Check schema for input validation https://github.com/keleshev/schema
# Check fire for exposing objects to CLI https://github.com/google/python-fire
# Check docopt for args parsing https://github.com/docopt/docopt

import json
from pathlib import Path
from docopt import docopt

# import fire


def validate_context(default, user):
    """
    Validates and merges defaults and user configurations.
    Received configeration dictionaries
    """
    # Ideas for validation:
    # number of base modes < number of snapshots
    # number of base mode > number of requested modes

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
        config_fname="configuration.json"

        try:
            context_user = json.loads((root_path / config_fname).read_text())[
                "config_data"
            ]
        except FileNotFoundError:
            #print("WARNING: No such configuration file: '{}'".format(config_fname))
            context_user = {}

        context_defaults = {
            # most frequently set
            "cases_test_dataset": [0],
            "rve_data_points": [300, 250, 200],
            "rve_data_points_range_list": [[100, 250, 50], [400, 600, 100]],
            "rve_data_points_rom": True,
            "rve_data_modes": [20, 30],
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
            # training files stuff
            "training_path": "sampling",
            "training_rve_materials_fname": "materials.json",
            "training_rve_model_fname": "model.mdpa",
            "case_path_pattern": "case_{}",
            "snapshots_fname": "snapshots.h5",
            "training_strain_fname": "_training_strain_set.dat",
            # offline files stuff
            "offline_path": "offline_data",
            "bases_fname_pattern": "bases_{}_{}m.npy",
            "local_bases_fname_pattern": "bases_inelastic_local_{}.npy",
            "local_sv_fname_pattern": "sv_inelastic_local_{}.dat",
            "roc_fname_pattern": "roc_{}ip",
            "rve_fname_pattern": "rve_{}m_{}ip.json",
            "correl_matrix_strain_pattern": "correlation_strain_{}.npy",
            "correl_matrix_damage_pattern": "correlation_r_value_{}.npy",
            # multiscale files stuff
            "multiscale_path": "multiscale_1ip",
            # other files stuff
        }
        config = validate_context(context_defaults, context_user)
        self.context = config
        # TODO:
        # load defaults
        # set default user location, overwrite with args
        # load user config
        # update default config with user config

        # file management
        self.root_path = root_path
        self.training_path = self.root_path / config["training_path"]
        self.offline_path = self.root_path / config["offline_path"]
        self.multiscale_path = self.root_path / config["multiscale_path"]

        #self.reuse_existing_files = config["reuse_existing_files"]
        #self.bases_fname = config["bases_fname_pattern"]
        #self.local_bases_fname = config["local_bases_fname_pattern"]
        #self.local_sv_fname = config["local_sv_fname_pattern"]

        # bases generation
        self.svd_cutoff = {}

        # self.energy_name = config["energy_name"]
        self.energy_elastic_modes = config["energy_elastic_modes"]
        self.energy_inelastic_modes = config["energy_inelastic_modes"]
        self.svd_cutoff[config["energy_name"]] = config["energy_svd_cutoff"]

        # self.strain_name = config["strain_name"]
        self.strain_elastic_modes = config["strain_elastic_modes"]
        self.strain_inelastic_modes = config["strain_inelastic_modes"]
        self.svd_cutoff[config["strain_name"]] = config["strain_svd_cutoff"]

        # self.rvalue_name = config["rvalue_name"]
        self.rvalue_elastic_modes = config["rvalue_elastic_modes"]
        self.rvalue_inelastic_modes = config["rvalue_inelastic_modes"]
        self.svd_cutoff[config["rvalue_name"]] = config["rvalue_svd_cutoff"]

        # points
        self.ip_subsets = [x for x in config["rve_data_points"]]
        for r in config["rve_data_points_range_list"]:  # unpack list of "ranges"
            for i in range(*r):
                self.ip_subsets.append(i)
        self.ip_subsets = list(set(self.ip_subsets))
        self.ip_subsets.sort()
        if config["rve_data_points_rom"]:
            self.ip_subsets.append("ROM")

        self.roc_fname_pattern = config["roc_fname_pattern"]

        # modes
        # self.reduced_nr_modes = config["rve_data_modes"]

        self.materials_fname = self.training_path / Path(
            config["training_rve_materials_fname"]
        )
        self.rve_fname_pattern = config["rve_fname_pattern"]

    def dump_config(self, fname):
        """
        Writes configuration file.
        """
        Path(fname).write_text(json.dumps({"config_data": self.context}, indent=2))

    def parse_training_strain_set(self):
        """
        Returns list of strain vectors used for trainig, read from file defined in configuration
        """
        fpath = self.training_path / self.context["training_strain_fname"]
        return fpath.read_text().splitlines()

    def case_name(self, c_id):
        """
        Returns case name with corresponging leading zeros
        e.g. if nr_cases:100, id:00..99, nr_id: 2 -> case_01..case_99
        """
        strain_set = self.parse_training_strain_set()
        nr_cases = len(strain_set)
        len_id = len(str(nr_cases - 1))  # size of the case number string
        case_id = "{:0{}d}".format(c_id, len_id)
        case_name = self.context["case_path_pattern"].format(case_id)
        return case_name

    def roc_fname(self, points):
        """
        docstrings here
        """
        return self.roc_fname_pattern.format(points)

    def rve_fname(self, modes, points):
        """
        docstrings here
        """
        return self.rve_fname_pattern.format(modes, points)

    def skip_calculation(self, fname):
        """ 
        Generates a list of files following filename pattern.
        Length of list is used as flag (False if empty, True otherwise)
        """
        fpath = Path.cwd() / fname  # converts filename to absolute Path
        flag_exists = len([f for f in fpath.parent.glob(fpath.name)])
        flag_reuse = self.context["reuse_existing_files"]
        return flag_exists and flag_reuse

    def get_bases_fname(self, field):
        """
        docstrings here
        """
        filename = self.context["bases_fname_pattern"].format(field, "*")
        fpath = self.offline_path / filename
        files = [f for f in fpath.parent.glob(fpath.name)]
        if len(files) == 0:
            return None
        if len(files) > 1:
            print(
                "Warning: More than one {} bases file detected. "
                "Picking first in the list: {}".format(field, files[0].name)
            )
        return files[0]


#####################################################################
# main
#####################################################################

if __name__ == "__main__":
    C = Common(root_path=Path("."))
    print("Test:")
    print(C.roc_fname("1"))
    print(C.roc_fname("100"))
    print(C.roc_fname("1000"))
    print(C.roc_fname(1000))
    print(C.roc_fname("ROM"))
    print(C.rve_fname(20, "1"))
    print(C.rve_fname(20, "100"))
    print(C.rve_fname(200, "1000"))
    print(C.rve_fname(200, 1000))
    print(C.rve_fname("2000", "ROM"))
    print(C.ip_subsets)
