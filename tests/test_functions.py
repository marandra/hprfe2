import pytest
from pathlib import Path
import numpy
from hprfe2.bases import Bases


b = Bases()
b_path = Path(b.config["local_bases_fname_pattern"].format("TEST"))
sv_path = Path(b.config["local_sv_fname_pattern"].format("TEST"))


def test_read_svd_01():
    """Test case when no cases are passed."""
    array = b.read_local_svd([], "TEST", 0.0001)
    numpy.testing.assert_allclose(array, numpy.empty([0, 0]))


def test_read_svd_02():
    """Test case when only one case is passed."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))
    out_array = numpy.array([[0.22, 0.12], [0.42, 0.22], [0.62, 0.32], [0.82, 0.42]])

    # Actual test
    array = b.read_local_svd(["."], "TEST", 0.0001)
    numpy.testing.assert_allclose(array, out_array)

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_read_svd_03():
    """Test case when more than one case are passed."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))
    out_array = numpy.array(
        [
            [0.22, 0.12, 0.22, 0.12],
            [0.42, 0.22, 0.42, 0.22],
            [0.62, 0.32, 0.62, 0.32],
            [0.82, 0.42, 0.82, 0.42],
        ]
    )
    # Actual test
    array = b.read_local_svd([".", "."], "TEST", 0.0001)
    numpy.testing.assert_allclose(array, out_array)

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_read_svd_04():
    """Test cutoff svd filter (only one mode passes)."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))
    out_array = numpy.array([[0.22], [0.42], [0.62], [0.82]])

    # Actual test
    array = b.read_local_svd(["."], "TEST", 1.5)
    numpy.testing.assert_allclose(array, out_array)

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_read_svd_05():
    """Test cutoff svd filter (no mode passes)."""

    # Generate test data
    in_array = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    numpy.save(b_path, in_array)
    numpy.savetxt(sv_path, numpy.array([2, 1]))

    # Actual test
    array = b.read_local_svd(["."], "TEST", 2.5)
    numpy.testing.assert_allclose(array, numpy.empty([4, 0]))

    # Remove test data
    b_path.unlink()
    sv_path.unlink()


def test_remove_elastic_modes_01():
    """Test removing elastic subspace from vectors, regular case"""

    # Generate test data
    in_vectors = numpy.array(
        [
            [0.11, 0.12, 0.13, 0.14],
            [0.21, 0.22, 0.23, 0.24],
            [0.31, 0.32, 0.33, 0.34],
            [0.41, 0.42, 0.43, 0.44],
        ]
    )
    in_subspace = numpy.array([[0.11, 0.12], [0.21, 0.22], [0.31, 0.32], [0.41, 0.42]])
    out_array = numpy.array(
        [
            [0.03506, 0.04262, 0.05018, 0.05774],
            [0.06994, 0.07538, 0.08082, 0.08626],
            [0.10482, 0.10814, 0.11146, 0.11478],
            [0.1397, 0.1409, 0.1421, 0.1433],
        ]
    )

    # Actual test
    array = b.remove_elastic_modes(in_vectors, in_subspace)
    numpy.testing.assert_allclose(array, out_array)


def test_remove_elastic_modes_02():
    """Test removing elastic subspace from vectors, null subspace"""

    # Generate test data
    in_vectors = numpy.array(
        [
            [0.11, 0.12, 0.13, 0.14],
            [0.21, 0.22, 0.23, 0.24],
            [0.31, 0.32, 0.33, 0.34],
            [0.41, 0.42, 0.43, 0.44],
        ]
    )
    in_subspace = numpy.zeros([4, 1])

    # Actual test
    array = b.remove_elastic_modes(in_vectors, in_subspace)
    numpy.testing.assert_allclose(array, in_vectors)


def test_remove_elastic_modes_03():
    """Test removing elastic subspace from vectors, no input vectors"""

    # Generate test data
    in_vectors = numpy.empty([4, 0])
    in_subspace = numpy.random.random([4, 1])

    # Actual test
    array = b.remove_elastic_modes(in_vectors, in_subspace)
    numpy.testing.assert_allclose(array, numpy.empty([4, 0]))


def test_compute_svd_01():
    """Test SVD with nr_modes=0"""

    # Generate test data
    in_vectors = numpy.array(
        [
            [0.11, 0.12, 0.13, 0.14],
            [0.21, 0.22, 0.23, 0.24],
            [0.31, 0.32, 0.33, 0.34],
            [0.41, 0.42, 0.43, 0.44],
        ]
    )

    # Actual test
    array = b.compute_svd(in_vectors, 0)
    numpy.testing.assert_allclose(array, numpy.empty([4, 0]))


def test_compute_svd_02():
    """Test standard full SVD with no vectors=0"""

    # Actual test
    array = b.compute_svd(numpy.ones([4, 0]), -1)
    numpy.testing.assert_allclose(array, numpy.ones([4, 0]))
