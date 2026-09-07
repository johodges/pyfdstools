"""Tests for reading FDS boundary (.bf) output."""

import numpy as np
import pytest

import pyfdstools as fds


def test_readBoundaryQuantities(case001_zip):
    quantities = fds.readBoundaryQuantities(case001_zip, 'case001')
    assert len(quantities) > 0
    assert 'WALL TEMPERATURE' in quantities


def test_queryBndf(case001_zip, data_dir):
    import os
    fdsFile = os.path.join(data_dir, 'case001.fds')
    datas, times = fds.queryBndf(
        case001_zip, 'case001', fdsFile, ['WALL TEMPERATURE'], -2, 4.4)
    assert 'WALL TEMPERATURE' in datas
    entry = datas['WALL TEMPERATURE']
    assert entry['DATA'].shape[:2] == entry['X'].shape == entry['Z'].shape
    assert entry['DATA'].shape[2] == len(times)
    assert np.all(np.diff(times) > 0)
    assert entry['UNITS'] == 'C'
    finite = entry['DATA'][np.isfinite(entry['DATA'])]
    assert finite.size > 0
    assert finite.min() > -50


def test_importBoundaryFile(case001_zip):
    bndfFiles = fds.getFileList(case001_zip, 'case001', 'bf')
    smvFile = fds.getFileList(case001_zip, 'case001', 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    times, patches, units = fds.importBoundaryFile(
        bndfFiles[0], smvFile=smvFile, gridNum=0, smvData=smvData)
    assert len(patches) > 0
    assert len(times) > 0
    for patch in patches:
        assert patch.data.shape[2] == len(times)


def test_readBoundaryHeader(case001_zip):
    bndfFile = fds.getFileList(case001_zip, 'case001', 'bf')[0]
    quantity, shortName, units, npatch = fds.readBoundaryHeader(bndfFile)
    assert isinstance(quantity, str) and len(quantity) > 0
    assert isinstance(units, str)
    assert npatch > 0


def test_fdspatch_buildSpace_shapes():
    """buildSpace must return coordinate grids matching the data grid.

    np.meshgrid defaults to 'xy' indexing, which transposes the result.
    That only goes unnoticed when the two patch dimensions are equal.
    """

    NX, NY, NT = 4, 7, 3
    lims = [1.0, 1.0, 0.0, 2.0, 0.0, 3.0]     # x is constant
    patch = fds.fdspatch(NX, NY, NT, lims, 1)
    patch.buildSpace()
    assert patch.x.shape == patch.y.shape == patch.z.shape == (NX+1, NY+1)
    # The constant axis holds the plane coordinate everywhere.
    assert np.allclose(patch.x, lims[0])
    # The in-plane axes vary along their own dimension only.
    assert np.allclose(np.diff(patch.y, axis=1), 0)
    assert np.allclose(np.diff(patch.z, axis=0), 0)


def test_fdspatch_average():
    patch = fds.fdspatch(3, 3, 4, [0.0, 0.0, 0.0, 1.0, 0.0, 1.0], 1)
    for i in range(0, 4):
        patch.append(np.full((3, 3), float(i)), i)
    assert np.allclose(patch.average([0, 1, 2, 3]), 1.5)
    assert np.allclose(patch.average([1, 2]), 1.5)


def test_getPatchLimsFromGrid_rejects_off_grid_bounds(case001_zip):
    """An unmatched patch bound must raise, not assert and return zeros."""
    smvFile = fds.getFileList(case001_zip, 'case001', 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    grid = smvData['grids'][0]
    from pyfdstools.extractBoundaryData import getPatchLimsFromGrid
    with pytest.raises(ValueError, match='do not lie on the mesh grid'):
        getPatchLimsFromGrid([1e9, 1e9, 1e9, 1e9, 1e9, 1e9], grid)


def test_buildAbsPatch_rejects_bad_axis(case001_zip):
    bndfFiles = fds.getFileList(case001_zip, 'case001', 'bf')
    smvFile = fds.getFileList(case001_zip, 'case001', 'smv')[0]
    smvData = fds.parseSMVFile(smvFile)
    _, patches, _ = fds.importBoundaryFile(
        bndfFiles[0], smvFile=smvFile, gridNum=0, smvData=smvData)
    with pytest.raises(ValueError, match='axis must be one of'):
        fds.buildAbsPatch(patches, 0, 1, 0, 1, 0, 1, 0.1, 0.1, 7)
