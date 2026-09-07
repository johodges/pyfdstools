"""Tests for reading FDS slice (.sf) output."""

import os

import numpy as np
import pytest

import pyfdstools as fds


# ---------------------------------------------------------------------
# Slice discovery
# ---------------------------------------------------------------------

def test_readSLCFquantities(case001_zip):
    quantities, files, dims, meshes, centers, units = \
        fds.readSLCFquantities('case001', case001_zip)
    assert len(quantities) == len(files) == len(dims) == len(units)
    assert 'TEMPERATURE' in quantities
    assert all(len(d) == 6 for d in dims)


def test_readSLCFquantities_returns_six_values(case001_zip):
    """paraview and the examples unpack six values from this call."""
    result = fds.readSLCFquantities('case001', case001_zip)
    assert len(result) == 6


def test_buildQTYstring(case001_zip):
    suffix = fds.buildQTYstring('case001', case001_zip, 'TEMPERATURE')
    assert isinstance(suffix, str)
    assert len(suffix) > 0


def test_buildQTYstring_unknown_quantity_raises(case001_zip):
    with pytest.raises(ValueError, match='not found'):
        fds.buildQTYstring('case001', case001_zip, 'NOT A REAL QUANTITY')


# ---------------------------------------------------------------------
# Reading a single slice file
# ---------------------------------------------------------------------

def slice_file(workingDir, chid):
    return sorted(fds.getFileList(workingDir, chid, 'sf'))[0]


def test_readSingleSlcfFile_all_times(case001_zip):
    sf = slice_file(case001_zip, 'case001')
    lims, data, times = fds.readSingleSlcfFile(sf)
    assert len(lims) == 6
    assert data.ndim == 4
    assert data.shape[3] == len(times)
    assert np.all(np.diff(times) > 0)
    assert np.isfinite(data).all()


def test_readSingleSlcfFile_single_time(case001_zip):
    sf = slice_file(case001_zip, 'case001')
    _, allData, times = fds.readSingleSlcfFile(sf)
    target = float(times[len(times)//2])
    _, oneFrame, oneTime = fds.readSingleSlcfFile(sf, time=target)
    assert oneFrame.shape[3] == 1
    assert np.isclose(oneTime[0], target)
    expected = allData[:, :, :, int(np.argmin(abs(np.asarray(times)-target)))]
    assert np.allclose(oneFrame[:, :, :, 0], expected)


def test_readSingleSlcfFile_time_average_matches_mean(case001_zip):
    """The averaged frame must be the mean of the frames in the window.

    Releases before v0.0.24 summed N frames and divided by N-1, which
    biased every averaged slice.
    """

    sf = slice_file(case001_zip, 'case001')
    _, allData, times = fds.readSingleSlcfFile(sf)
    times = np.asarray(times)

    query, window = 30.0, 20.0
    _, averaged, avgTime = fds.readSingleSlcfFile(sf, time=query, dt=window)

    inds = np.where((times >= query - window/2)
                    & (times <= query + window/2))[0]
    assert inds.size > 1, "test needs several frames inside the window"
    expected = allData[:, :, :, inds].mean(axis=3)

    assert np.allclose(averaged[:, :, :, 0], expected, rtol=1e-5, atol=1e-4)
    assert np.isclose(avgTime[0], (times[inds[0]] + times[inds[-1]])/2)


def test_readSingleSlcfFile_matches_between_zip_and_directory(
        case001_zip, case001_dir):
    """Reading from an archive and from a directory must agree.

    The two take different code paths: seek-based reads for a real file,
    stream reads for an archive member.
    """

    zipFile = slice_file(case001_zip, 'case001')
    dirFile = slice_file(case001_dir, 'case001')
    assert os.path.basename(zipFile) == os.path.basename(dirFile)

    limsZip, dataZip, timesZip = fds.readSingleSlcfFile(zipFile)
    limsDir, dataDir, timesDir = fds.readSingleSlcfFile(dirFile)

    assert limsZip == limsDir
    assert np.allclose(np.asarray(timesZip), np.asarray(timesDir))
    assert np.allclose(dataZip, dataDir)


def test_readSLCFtimes_matches_frame_count(case001_zip):
    sf = slice_file(case001_zip, 'case001')
    times = fds.readSLCFtimes(sf)
    _, data, _ = fds.readSingleSlcfFile(sf)
    assert len(times) == data.shape[3]


def test_readSLCFtimes_caches_to_file(case001_dir, outdir):
    sf = slice_file(case001_dir, 'case001')
    cache = os.path.join(outdir, 'times.csv')
    first = fds.readSLCFtimes(sf, timesFile=cache)
    assert os.path.exists(cache)
    second = fds.readSLCFtimes(sf, timesFile=cache)
    assert np.allclose(first, second)


# ---------------------------------------------------------------------
# 2-D slice queries
# ---------------------------------------------------------------------

def test_query2dAxisValue(case001_zip):
    data, units = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    assert data is not None
    assert units == 'C'
    assert data['datas'].shape[:2] == data['x'].shape == data['z'].shape
    assert data['datas'].shape[2] == len(data['times'])
    assert np.isfinite(data['datas']).all()
    # An FDS temperature slice of a compartment fire is physical.
    assert data['datas'].min() > -50
    assert data['datas'].max() < 2000


def test_query2dAxisValue_all_times(case001_zip):
    data, units = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55)
    assert data['datas'].shape[2] > 1
    assert np.all(np.diff(data['times']) > 0)


def test_query2dAxisValue_missing_plane_returns_none(case001_zip, capsys):
    data, units = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 3, 999.0)
    assert data is None and units is None
    printed = capsys.readouterr().out
    assert 'Available slices' in printed


def test_query2dAxisValue_coordinates_are_monotonic(case001_zip):
    data, _ = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    assert np.all(np.diff(data['x'][:, 0]) > 0)
    assert np.all(np.diff(data['z'][0, :]) > 0)


def test_read2dSliceFile_running_average(case001_zip):
    """A running average must equal the direct per-frame mean.

    The prefix-sum implementation is an optimisation of exactly that
    calculation, so it has to reproduce it.
    """

    sf = slice_file(case001_zip, 'case001')
    x, z, plain, times, coords = fds.read2dSliceFile(sf, 'case001')
    times = np.asarray(times)
    window = 20.0
    _, _, averaged, _, _ = fds.read2dSliceFile(sf, 'case001', dt=window)

    assert averaged.shape == plain.shape
    for i in (0, len(times)//2, len(times)-1):
        inds = np.where((times >= times[i] - window/2)
                        & (times <= times[i] + window/2))[0]
        expected = plain[:, :, inds].mean(axis=2)
        assert np.allclose(averaged[:, :, i], expected, rtol=1e-5, atol=1e-4)


def test_read2dSliceFile_coords(case001_zip):
    sf = slice_file(case001_zip, 'case001')
    x, z, data, times, coords = fds.read2dSliceFile(sf, 'case001')
    assert len(coords) == 6
    assert x.shape == z.shape == data.shape[:2]
    # A 2-D slice is flat along exactly one axis.
    flat = [coords[0] == coords[1], coords[2] == coords[3],
            coords[4] == coords[5]]
    assert sum(flat) == 1


# ---------------------------------------------------------------------
# 3-D slice queries
# ---------------------------------------------------------------------

def test_readSLCF3Ddata(case001_zip):
    grid, data, times, units = fds.readSLCF3Ddata(
        'case001', case001_zip, 'TEMPERATURE')
    assert grid.shape[3] == 3
    assert data.shape[:3] == grid.shape[:3]
    assert data.shape[3] == len(times)
    assert units == 'C'


def test_readSLCF3Ddata_unknown_quantity(case001_zip):
    """A quantity with no 3-D slice must report itself, not crash."""
    grid, data, times, units = fds.readSLCF3Ddata(
        'case001', case001_zip, 'NOT A REAL QUANTITY')
    assert grid is False and data is False


def test_readSLCF3Ddata_plane_extraction(case001_zip):
    grid, data, times, units = fds.readSLCF3Ddata(
        'case001', case001_zip, 'TEMPERATURE', axis=3, value=1.0)
    assert grid.ndim == 3 and grid.shape[2] == 2
    assert data.ndim == 3


def test_findSliceLocation(case001_zip):
    grid, data, times, units = fds.readSLCF3Ddata(
        'case001', case001_zip, 'TEMPERATURE')
    x, z, slc = fds.findSliceLocation(grid, data[:, :, :, -1], 3, 1.0)
    assert x.shape == z.shape == slc.shape
