"""
BASES: pending description here.
"""
import logging
import os
import time
from pathlib import Path
# import multiprocessing
import numpy
import h5py
import sklearn.decomposition
from common import Common


logger = logging.getLogger(__name__)


class Bases(Common):
    """Functions for snapshots I/O and bases generation."""

    def __init__(self, **kargv):
        super().__init__(**kargv)

    def write_field_to_hdf5(self, filename, group, field, timestep, data):
        with h5py.File(filename, "a") as f:
            f.create_dataset("{}/{}/{}".format(group, field, timestep), data=data)

    def _get_shape_of_snapshots_in_case(self, spath, group, field):
        """
        Receives path and filename of snapshots file
        """
        rows = 0
        cols = 0
        with h5py.File(spath, "r") as f:
            try:
                d = f[group][field]
                for k, v in d.items():
                    rows = len(v)  # TODO not optimal. read only once
                    cols += 1
            except KeyError:
                # not counting missing datasets
                pass
        return rows, cols

    def _read_snapshots_in_case(self, spath, group, field):
        """
        Receives path and filename of snapshots file
        """
        rows, cols = self._get_shape_of_snapshots_in_case(spath, group, field)
        snapshots = numpy.empty([rows, cols])
        column = 0
        with h5py.File(spath, "r") as f:
            try:
                d = f[group][field]
                for k, v in d.items():
                    snapshots[:, column] = v
                    column += 1
            except KeyError:
                logger.debug(
                    "    skipping {}/{} of {} (dataset not present)".format(
                        group, field, spath.parent.name
                    )
                )
        return snapshots

    def read_snapshots(self, cases, group, field):
        fname = self.config["snapshots_fname"]
        paths = sorted([f / fname for f in cases])

        logger.debug("  - getting shape of snapshots to allocate array")
        rows = 0
        cols = 0
        for path in paths:
            r, c = self._get_shape_of_snapshots_in_case(path, group, field)
            cols += c
            rows = r

        logger.info("  - loading {} snapshots".format(cols))
        arrays = numpy.empty([rows, cols])
        batch_size = int(len(paths) / 10 + 0.5)
        counter = 1
        column = 0
        for path in paths:
            array = self._read_snapshots_in_case(path, group, field)
            if numpy.shape(array)[1] == 0:  # missing dataset
                continue
            arrays[:, column : column + numpy.shape(array)[1]] = array
            column += numpy.shape(array)[1]
            #
            if not counter % batch_size:
                logger.info(
                    "    {}/{} trajectories processed".format(counter, len(paths))
                )
            counter += 1
        return arrays

    def read_local_svd(self, cases, field, cutoff_tol):
        """Return array combining bases of training cases.

        Return array combining bases precomputed (svd) on each training case.
        Modes to be included must have a corresponding singular value greater than cutoff value.

        Arguments:
            cases {list of str or Path} -- list with locations where to look for bases
            field {str} -- field name
            cutoff_tol {float} -- singular value below which the corresponfing mode from base wnt be loaded

        Returns:
            numpy.arrary -- combined array with bases from the training cases that fit the cutoff criteria
        """
        # TODO: Add test for feature ignore missing local SVD.
        b_fname = self.config["local_bases_fname_pattern"].format(field)
        sv_fname = self.config["local_sv_fname_pattern"].format(field)
        paths = sorted([Path(f) for f in cases])

        # First run: compute space of the final array
        logger.debug("  - getting shape of local bases to allocate array")
        rows = 0
        cols = 0
        for path in paths:
            if not (
                path / b_fname
            ).exists():  # in case there are no inelastic snapshots, local bases are not generated and no local bases file present
                continue
            a = numpy.load(str(path / b_fname), mmap_mode="r")
            if numpy.shape(a)[1] == 0:  # missing dataset
                continue
            sv = numpy.loadtxt(path / sv_fname)
            idx = numpy.where(sv > cutoff_tol)[0]
            c = len(idx)
            cols += c
            rows = numpy.shape(a)[0]

        # Second run: load bases in array
        logger.info("  - loading {} inelastic modes".format(cols))
        arrays = numpy.empty([rows, cols])
        batch_size = int(len(paths) / 10 + 1.0)
        counter = 1
        column = 0
        for path in paths:
            if not (
                path / b_fname
            ).exists():  # in case there are no inelastic snapshots, local bases are not generated and no local bases file present
                continue
            a = numpy.load(str(path / b_fname), mmap_mode="r")
            if numpy.shape(a)[1] == 0:  # missing dataset
                continue
            sv = numpy.loadtxt(path / sv_fname)
            # This is a workaround. When size==1 we are getting sv=12.3,
            # instead of sv=[12.3] (that we expect)
            if sv.size == 1:
                sv = sv.reshape([1])
            idx = numpy.where(sv > cutoff_tol)[0]
            c = len(idx)
            arrays[:, column : column + c] = a[:, idx] * sv[idx]
            column += c

            if not counter % batch_size:
                logger.info(
                    "    {}/{} trajectories processed".format(counter, len(paths))
                )
            counter += 1

        return arrays

    def remove_elastic_modes(self, X, Ue):
        """Remove components in the base U from the vectors X.

        Function tested.

        Arguments:
            X {numpy.array} -- vectors to remove a subespace from
            Ue {numpy.array} -- subspace to remove

        Returns:
            numpy.array -- X without its proxection in Ue
        """
        logger.info("Removing elastic componennt")
        t0 = time.time()
        for i in range(numpy.shape(X)[1]):
            projection = X[:, i] @ Ue
            X[:, i] -= numpy.sum(projection * Ue, axis=1)
        logger.debug("    elapsed time: {:.1f}s".format(time.time() - t0))
        return X

    def compute_svd(self, X, nr_modes):
        """Compute (truncated) SVD decomposition of X vectors.

        Function tested.

        Arguments:
            X {numpy.arrar} -- Vectors to decompose
            nr_modes {int} -- If >= 0, then a partial (truncated) decomposition os performed using randomized algorithm.
                              If > 0, the a full decompisition is performed.

        Returns:
            numpy.array -- Bases. Writes file with singular values.
        """
        t0 = time.time()
        if nr_modes > -1:
            logger.info("- Computing SVD using RANDOMIZED algorithm")
            svd = sklearn.decomposition.TruncatedSVD(
                n_components=nr_modes, algorithm="randomized"
            )
            svd.fit(X.T)
            U = svd.components_.T
            S = svd.singular_values_.T
        else:
            logger.info("- Computing SVD using STANDARD algorithm")
            [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]

        logger.info("    - SVD time: {:.1f}s".format(time.time() - t0))
        logger.info("    - singular value of selected modes:")
        logger.info("      {}".format(S[:nr_modes]))
        logger.info("      validation: following singular values (excluded):")
        logger.info("      {}".format(S[nr_modes : nr_modes + 4]))
        logger.info("    - nr and size of modes: {}, {}".format(U.shape[1], U.shape[0]))
        logger.info("")
        numpy.savetxt(self.bases_path / "singular_values.dat", S)
        return U

    def create_bases(
        self,
        field_name,
        nr_elastic_modes,
        nr_inelastic_modes,
        cases_path,
        bases_fname,
        cutoff_tol,
    ):
        if self.skip_calculation(self.bases_path / bases_fname.format(field_name, "*")):
            logger.info(
                "File {} exists. Skipping calculation".format(
                    bases_fname.format(field_name, "*")
                )
            )
            return
        logger.info("Generating {} bases".format(field_name))

        t0 = time.time()
        # Snapshots splitted in elastic and inelastic groups
        if nr_elastic_modes > 0:
            logger.info("- Processing ELASTIC snapshots")
            X = self.read_snapshots(cases_path, "ELASTIC", field_name)
            Ue = self.compute_svd(X, nr_elastic_modes)
            os.rename(
                self.bases_path / "singular_values.dat",
                self.bases_path / "sv_{}_elastic.dat".format(field_name),
            )

            logger.info("- Processing INELASTIC modes")
            X = self.read_local_svd(cases_path, field_name, cutoff_tol)
            X = self.remove_elastic_modes(X, Ue)
            Ui = self.compute_svd(X, nr_inelastic_modes)
            os.rename(
                self.bases_path / "singular_values.dat",
                self.bases_path / "sv_{}_inelastic.dat".format(field_name),
            )

            U = numpy.hstack([Ue, Ui])

        # No splitting of elastic and inelastic snapshots
        else:
            logger.info(
                "Nr of elastic modes set to zero -> "
                "Not discriminating elastic/inelastic snapshots"
            )
            X = self.read_snapshots(cases_path, "INELASTIC", field_name)
            U = self.compute_svd(X, nr_inelastic_modes)
            os.rename(
                self.bases_path / "singular_values.dat",
                self.bases_path / "sv_{}.dat".format(field_name),
            )

        numpy.save(
            self.bases_path / bases_fname.format(field_name, numpy.shape(U)[1]), U
        )
        logger.info("  Elapsed time: {:.1f}s".format(time.time() - t0))
        logger.info("")

    def generate_local_bases(self, case, field, ss_fname, lb_fname, sv_fname):
        base = case / lb_fname
        logger.debug("   - missing {} {}".format(base.parent.name, base.name))
        X = self._read_snapshots_in_case(case / ss_fname, "INELASTIC", field)
        [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]
        numpy.save(base, U)
        path = case / sv_fname
        numpy.savetxt(path, S)

    def generate_missing_local_bases(self, field, threads=1):
        logger.info("Looking for missing local bases {}".format(field))
        cases_path = self.training_path.glob(
            self.config["case_path_pattern"].format("*")
        )
        lb_fname = self.config["local_bases_fname_pattern"].format(field)
        sv_fname = self.config["local_sv_fname_pattern"].format(field)
        ss_fname = self.config["snapshots_fname"]
        missing = []
        for case in cases_path:
            bases = case / lb_fname
            if self.skip_calculation(bases):
                continue
            missing.append(case)
        if not missing:
            return
        # There are missing bases files. Let's generate them.
        for case in missing:
            try:
                self.generate_local_bases(case, field, ss_fname, lb_fname, sv_fname)
            except:
                logger.info("  - skipping case: not inelastic snapshots present")
                continue

        # Testing: version with Pool
        # with multiprocessing.Pool(processes=threads) as pool:
        #    logger.debug("   - generating bases")
        #    logger.debug("   - multiprocessing {} threads".format(threads))
        #    pool.starmap(generate_local_bases, zip(missing, [field] * len(missing)))

        # Testing: version with Process
        # processes = []
        # for case in missing:
        #    semaphore.acquire()
        #    p = multiprocessing.Process(target=generate_local_bases, args=(case, field, semaphore))
        #    processes.append(p)
        #    p.start()
        # for p in processes:
        #    p.join()


def run(common):
    """Creates file structure from the computation of bases"""

    logger.info("Beginning bases calculation -----------------------")

    #
    # removing cases from training dataset
    # TODO: add TRAINING set and TEST set as members of Common
    #
    training_set = []
    for c in common.training_path.glob(B.config["case_path_pattern"].format("*")):
        c_id = int(c.name.split("_")[1])
        if c_id in common.config["cases_test_dataset"]:
            logger.info("Removing case {} from training dataset".format(c.name))
            continue
        training_set.append(c)

    #
    # generate missing local bases
    #
    B.generate_missing_local_bases(
        common.config["energy_name"],
    )
    B.generate_missing_local_bases(
        common.config["strain_name"],
    )
    B.generate_missing_local_bases(
        common.config["rvalue_name"],
    )

    #
    # compute bases
    #
    B.create_bases(
        common.config["energy_name"],
        common.config["energy_elastic_modes"],
        common.config["energy_inelastic_modes"],
        training_set,
        common.config["bases_fname_pattern"],
        common.svd_cutoff[common.config["energy_name"]],
    )
    B.create_bases(
        common.config["strain_name"],
        common.config["strain_elastic_modes"],
        common.config["strain_inelastic_modes"],
        training_set,
        common.config["bases_fname_pattern"],
        common.svd_cutoff[common.config["strain_name"]],
    )
    B.create_bases(
        common.config["rvalue_name"],
        common.config["rvalue_elastic_modes"],
        common.config["rvalue_inelastic_modes"],
        training_set,
        common.config["bases_fname_pattern"],
        common.svd_cutoff[common.config["rvalue_name"]],
    )
    logger.info("Finished -----------------------------------------")

    return
