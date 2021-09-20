"""
BASES: Functions for snapshots I/O and bases generation.
"""
import logging
import os
import time
from pathlib import Path

import h5py
# import multiprocessing
import numpy
import sklearn.decomposition

logger = logging.getLogger(__name__)


def write_field_to_hdf5(filename, group, field, timestep, data):
    with h5py.File(filename, "a") as f:
        f.create_dataset("{}/{}/{}".format(group, field, timestep), data=data)


def _get_shape_of_snapshots_in_case(spath, group, field):
    """
    Receives path and filename of snapshots file
    """
    rows = 0
    cols = 0
    with h5py.File(spath, "r") as f:
        try:
            d = f[group][field]
            for _, v in d.items():
                rows = len(v)  # TODO not optimal. read only once
                cols += 1
        except KeyError:
            pass
    return rows, cols


def _read_snapshots_in_case(spath, group, field):
    """
    Receives path and filename of snapshots file
    """
    rows, cols = _get_shape_of_snapshots_in_case(spath, group, field)
    snapshots = numpy.empty([rows, cols])
    column = 0
    with h5py.File(spath, "r") as f:
        try:
            d = f[group][field]
            for _, v in d.items():
                snapshots[:, column] = v
                column += 1
        except KeyError:
            logger.debug(f"    skipping {group}/{field} of {spath.parent.name} (dataset not present)")
    return snapshots


def read_snapshots(common, cases, group, field):
    fname = common.config["snapshots_fname"]
    paths = sorted([f / fname for f in cases])

    logger.debug("  - getting shape of snapshots to allocate array")
    rows = 0
    cols = 0
    for path in paths:
        r, c = _get_shape_of_snapshots_in_case(path, group, field)
        cols += c
        rows = max(r, rows)

    logger.info("  - loading {} snapshots".format(cols))
    arrays = numpy.empty([rows, cols])
    batch_size = int(len(paths) / 10 + 0.5)
    counter = 1
    column = 0
    for path in paths:
        array = _read_snapshots_in_case(path, group, field)
        if numpy.shape(array)[1] == 0:  # missing dataset
            continue
        arrays[:, column : column + numpy.shape(array)[1]] = array
        column += numpy.shape(array)[1]
        #
        if not counter % batch_size:
            logger.info("    {}/{} trajectories processed".format(counter, len(paths)))
        counter += 1
    return arrays


def read_local_svd(common, cases, field, cutoff_tol):
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
    b_fname = common.config["local_bases_fname_pattern"].format(field)
    sv_fname = common.config["local_sv_fname_pattern"].format(field)
    paths = sorted([Path(f) for f in cases])

    # First run: compute space of the final array
    logger.debug("  - getting shape of local bases to allocate array")
    rows = 0
    cols = 0
    for path in paths:
        # in case there are no inelastic snapshots, local bases are not
        # generated and no local bases file present
        if not (path / b_fname).exists():
            continue
        logger.debug(f"        {str(path / b_fname)}")
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
            logger.info("    {}/{} trajectories processed".format(counter, len(paths)))
        counter += 1

    return arrays


def remove_elastic_modes(X, Ue):
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


def compute_svd(common, X, nr_modes):
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
    numpy.savetxt(common.bases_path / "singular_values.dat", S)
    return U


def create_bases(
    common,
    field,
    cases_path,
):
    logger.info(f"Generating {field} bases")

    t0 = time.time()

    nr_elastic_modes = common.config[field]["nr_mode_elastic"]
    nr_inelastic_modes = common.config[field]["nr_mode_inelastic"]
    cutoff_tol = common.config[field]["svd_cutoff"]

    # Snapshots splitted in elastic and inelastic groups
    if nr_elastic_modes > 0:
        logger.info("- Processing ELASTIC snapshots")
        X = read_snapshots(common, cases_path, "ELASTIC", field)
        Ue = compute_svd(common, X, nr_elastic_modes)
        os.rename(
            common.bases_path / "singular_values.dat",
            common.bases_path / f"sv_{field}_elastic.dat",
        )

        logger.info("- Processing INELASTIC modes")
        X = read_local_svd(common, cases_path, field, cutoff_tol)
        X = remove_elastic_modes(X, Ue)
        Ui = compute_svd(common, X, nr_inelastic_modes)
        os.rename(
            common.bases_path / "singular_values.dat",
            common.bases_path / f"sv_{field}_inelastic.dat",
        )

        U = numpy.hstack([Ue, Ui])

    # No splitting of elastic and inelastic snapshots
    else:
        logger.info(
            "Nr of elastic modes set to zero -> "
            "Not discriminating elastic/inelastic snapshots"
        )
        X = read_snapshots(cases_path, "INELASTIC", field)
        U = compute_svd(X, nr_inelastic_modes)
        os.rename(
            common.bases_path / "singular_values.dat",
            common.bases_path / f"sv_{field}.dat",
        )

    logger.info("  Elapsed time: {:.1f}s".format(time.time() - t0))
    logger.info("")
    return U


def generate_local_bases(case, field, ss_fname, lb_fname, sv_fname):
    base = case / lb_fname
    logger.debug("   - missing {} {}".format(base.parent.name, base.name))
    X = _read_snapshots_in_case(case / ss_fname, "INELASTIC", field)
    [U, S] = numpy.linalg.svd(X, full_matrices=False)[:2]
    numpy.save(base, U)
    path = case / sv_fname
    numpy.savetxt(path, S)


def generate_missing_local_bases(common, field, threads=1):
    logger.info("Looking for missing local bases {}".format(field))
    cases_path = common.training_path.glob(
        common.config["sampling_case_path_pattern"].format("*")
    )
    lb_fname = common.config["local_bases_fname_pattern"].format(field)
    sv_fname = common.config["local_sv_fname_pattern"].format(field)
    ss_fname = common.config["snapshots_fname"]
    missing = []
    for case in cases_path:
        bases = case / lb_fname
        if common.skip_calculation(bases):
            continue
        missing.append(case)
    if not missing:
        return
    # There are missing bases files. Let's generate them.
    for case in missing:
        try:
            generate_local_bases(case, field, ss_fname, lb_fname, sv_fname)
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

    training_set = []
    for c in common.training_path.glob(common.config["sampling_case_path_pattern"].format("*")):
        training_set.append(c)

    #
    # generate missing local bases
    #
    for field in ["ENERGY", "STRAIN", "RVALUE"]:
        generate_missing_local_bases(
            common,
            field,
        )

    #
    # compute bases
    #
    for field in ["ENERGY", "STRAIN", "RVALUE"]:
        if common.has_dataset("BASES", field):
            logger.info(f"Dataset BASES/{field} exits. Skipping.")
        else:
            U = create_bases(
                common,
                field,
                training_set,
            )
            common.set_dataset(U, "BASES", field)

    return
