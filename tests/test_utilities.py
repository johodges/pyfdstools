"""Tests for the shared utility routines."""

import os

import numpy as np
import pytest

import pyfdstools as fds


# ---------------------------------------------------------------------
# File discovery and IO helpers
# ---------------------------------------------------------------------

def test_getFileList_zip_and_dir_agree(case001_zip, case001_dir):
    fromZip = fds.getFileList(case001_zip, 'case001', 'sf')
    fromDir = fds.getFileList(case001_dir, 'case001', 'sf')
    assert len(fromZip) == len(fromDir) > 0
    assert (sorted(os.path.basename(x) for x in fromZip)
            == sorted(os.path.basename(x) for x in fromDir))


def test_getFileListFromResultDir_is_an_alias(case001_zip):
    assert (sorted(fds.getFileListFromResultDir(case001_zip, 'case001', 'sf'))
            == sorted(fds.getFileList(case001_zip, 'case001', 'sf')))


def test_getFileList_unknown_extension(case001_zip):
    assert fds.getFileList(case001_zip, 'case001', 'nosuchext') == []


def test_zopen_reads_from_archive_and_directory(case001_zip, case001_dir):
    zipFile = fds.getFileList(case001_zip, 'case001', 'smv')[0]
    dirFile = fds.getFileList(case001_dir, 'case001', 'smv')[0]
    assert fds.zreadlines(zipFile) == fds.zreadlines(dirFile)


def test_getEndianness(case001_zip):
    assert fds.getEndianness(case001_zip, 'case001') in ('<', '>')


def test_getDatatypeByEndianness():
    assert fds.getDatatypeByEndianness(np.float32, '<').byteorder in ('<', '=')
    assert fds.getDatatypeByEndianness(np.float32, '>').byteorder == '>'
    with pytest.raises(ValueError):
        fds.getDatatypeByEndianness(np.float32, '?')


# ---------------------------------------------------------------------
# Grid helpers
# ---------------------------------------------------------------------

def test_getAbsoluteGrid_spans_every_mesh(case002_zip):
    smvFile = fds.getFileList(case002_zip, 'case002', 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    grids = {}
    for i, (trnx, trny, trnz) in enumerate(smvData['grids']):
        xGrid, yGrid, zGrid = np.meshgrid(
            trnx[:, 1], trny[:, 1], trnz[:, 1], indexing='ij')
        grids['mesh%d' % (i)] = {
            'xGrid': xGrid, 'yGrid': yGrid, 'zGrid': zGrid}
    absGrid = fds.getAbsoluteGrid(grids)
    assert absGrid.shape[3] == 3
    for key in grids:
        assert absGrid[:, :, :, 0].min() <= grids[key]['xGrid'].min()
        assert absGrid[:, :, :, 0].max() >= grids[key]['xGrid'].max()


def test_pointsFromXB():
    pts = fds.pointsFromXB([0, 1, 0, 2, 0, 3])
    assert np.asarray(pts).shape == (8, 3)
    assert np.asarray(pts).min() == 0
    assert np.asarray(pts).max() == 3


def test_pointsFromXB_extend():
    pts = np.asarray(fds.pointsFromXB([0, 1, 0, 1, 0, 1], extend=[0.5]*3))
    assert pts.min() == -0.5 and pts.max() == 1.5


# ---------------------------------------------------------------------
# Signal processing
# ---------------------------------------------------------------------

def test_kalmanFilter_smooths_noise():
    rng = np.random.default_rng(0)
    truth = 5.0
    noisy = truth + rng.normal(0, 0.5, 500)
    filtered = fds.kalmanFilter(noisy)
    assert filtered.shape == noisy.shape
    assert filtered[100:].std() < noisy.std()
    assert abs(filtered[-1] - truth) < 0.5


def test_timeAverage_of_a_constant_is_the_constant():
    times = np.linspace(0, 100, 201)
    data = np.full((3, 4, times.size), 7.0)
    averaged, outTimes = fds.timeAverage(data, times, 10.0)
    assert np.allclose(averaged, 7.0)
    assert averaged.shape[:2] == (3, 4)
    assert averaged.shape[2] == len(outTimes)


def test_timeAverage_reduces_oscillation():
    """A boxcar wider than the period must average the oscillation out.

    With smoothEnds=False the first and last half-window are copied
    through unaveraged by design, so only the interior is checked.
    """

    times = np.linspace(0, 100, 1001)
    window = 10.0
    signal = 10.0 + np.sin(2*np.pi*times)
    data = signal.reshape(1, 1, -1)
    averaged, outTimes = fds.timeAverage(data, times, window)

    halfWindow = int(round(window/(times[1] - times[0])))
    interior = averaged[0, 0, halfWindow:-halfWindow]
    assert interior.size > 0
    assert abs(interior.mean() - 10.0) < 0.05
    assert interior.std() < 0.1
    assert data[0, 0].std() > 0.5


def test_timeAverage_leaves_ends_unaveraged_by_default():
    """smoothEnds=False copies the half-window at each end through."""
    times = np.linspace(0, 100, 1001)
    signal = 10.0 + np.sin(2*np.pi*times)
    data = signal.reshape(1, 1, -1)
    plain, _ = fds.timeAverage(data, times, 10.0)
    smoothed, _ = fds.timeAverage(data, times, 10.0, smoothEnds=True)
    assert plain.shape == smoothed.shape
    assert not np.allclose(plain[0, 0, :10], smoothed[0, 0, :10])


def test_timeAverage_window_larger_than_series_is_a_noop():
    times = np.linspace(0, 10, 11)
    data = np.ones((1, 1, 11))
    averaged, outTimes = fds.timeAverage(data, times, 100.0)
    assert np.allclose(averaged, data)
    assert np.allclose(outTimes, times)


def test_timeAverage2_is_deprecated():
    times = np.linspace(0, 100, 201)
    data = np.full((2, 2, times.size), 3.0)
    with pytest.warns(DeprecationWarning):
        result = fds.timeAverage2(data, times, 10.0)
    assert result.shape == data.shape


# ---------------------------------------------------------------------
# Two-zone reduction
# ---------------------------------------------------------------------

def test_getTwoZone_uniform_profile():
    """A uniform profile has no layer interface to find."""
    z = np.linspace(0, 3, 50)
    val = np.full_like(z, 300.0)
    low, high, interface = fds.getTwoZone(z, val)
    assert np.isclose(low, 300.0)
    assert np.isclose(high, 300.0)
    assert np.isclose(interface, z.max())


def test_getTwoZone_hot_upper_layer():
    """A hot upper layer must be detected above a cool lower layer."""
    z = np.linspace(0, 3, 200)
    val = np.where(z > 2.0, 500.0, 25.0)
    low, high, interface = fds.getTwoZone(z, val)
    assert high > low
    assert 1.0 < interface < 3.0
    assert 20.0 < low < 120.0
    assert high > 300.0


def test_getTwoZone_runs_on_numpy_2():
    """np.trapz was removed in numpy 2.0; getTwoZone must still work."""
    z = np.linspace(0, 3, 50)
    val = 25.0 + 400.0*np.exp(-((z - 2.6)**2)/0.4)
    result = fds.getTwoZone(z, val)
    assert len(result) == 3
    assert all(np.isfinite(result))


# ---------------------------------------------------------------------
# Adiabatic surface temperature
# ---------------------------------------------------------------------

def test_astFromGhf_round_trips():
    """AST must satisfy the gauge energy balance it was solved from."""
    sigma = 5.67e-11
    h, e, Tgauge = 0.01, 0.9, 20.0
    for ghf in (1.0, 5.0, 20.0):
        ast = fds.astFromGhf(ghf, h, e, Tgauge=Tgauge)
        residual = (e*sigma*((ast + 273.15)**4 - (Tgauge + 273.15)**4)
                    + h*(ast - Tgauge) - ghf)
        assert abs(residual) < 1e-6*max(1.0, abs(ghf))


def test_astFromGhf_is_monotonic_in_flux():
    h, e = 0.01, 0.9
    fluxes = np.array([1.0, 2.0, 5.0, 10.0, 20.0])
    asts = np.array([fds.astFromGhf(f, h, e) for f in fluxes])
    assert np.all(np.diff(asts) > 0)


# ---------------------------------------------------------------------
# Colors
# ---------------------------------------------------------------------

def test_buildSMVcolormap():
    cmap = fds.buildSMVcolormap()
    assert cmap.N > 0
    rgba = cmap(0.5)
    assert len(rgba) == 4


def test_getPlotColors():
    colors = fds.getPlotColors(5)
    assert len(colors) == 5
    assert all(c.startswith('#') for c in colors)
