"""Tests for parsing FDS smokeview (.smv) files."""

import numpy as np
import pytest

import pyfdstools as fds


def parse(workingDir, chid):
    smvFile = fds.getFileList(workingDir, chid, 'smv')[0]
    return fds.parseSMVFile(smvFile)


def test_single_mesh_case(case001_zip):
    smv = parse(case001_zip, 'case001')
    assert len(smv['grids']) == 1
    gridTRNX, gridTRNY, gridTRNZ = smv['grids'][0]
    # Each TRN table has one row per node, so NX + 1 rows for NX cells.
    assert gridTRNX.shape[1] == 2
    assert gridTRNX.shape[0] > 1
    # Coordinates must increase monotonically along each axis.
    for table in (gridTRNX, gridTRNY, gridTRNZ):
        assert np.all(np.diff(table[:, 1]) > 0)


def test_multi_mesh_case(case002_zip):
    smv = parse(case002_zip, 'case002')
    assert len(smv['grids']) == 4
    assert len(smv['files']['SLICES']) > 0


def test_obstructions_have_finite_extent(case001_zip):
    smv = parse(case001_zip, 'case001')
    obsts = np.asarray(smv['obsts'])
    assert obsts.shape[0] > 0
    # Columns 13-18 are the bounding box in coordinates. parseOBST
    # extends thin obstructions by one cell, so every axis must span a
    # positive distance.
    for lo, hi in ((13, 14), (15, 16), (17, 18)):
        assert np.all(obsts[:, hi] > obsts[:, lo])


def test_stretched_mesh_parses(stretched_zip):
    """A stretched mesh inserts extra lines into the TRN blocks.

    Locating the coordinate tables at a fixed offset from the GRID line
    raised a ValueError on these cases.
    """

    smv = parse(stretched_zip, 'stretched_mesh_example')
    assert len(smv['grids']) == 14
    for gridTRNX, gridTRNY, gridTRNZ in smv['grids']:
        for table in (gridTRNX, gridTRNY, gridTRNZ):
            assert table.ndim == 2
            assert table.shape[1] == 2
            assert np.all(np.diff(table[:, 1]) > 0)

    # The z-axis of this case is genuinely stretched, so its cell sizes
    # must not all be equal.
    deltas = np.diff(smv['grids'][0][2][:, 1])
    assert deltas.max() - deltas.min() > 1e-6


def test_cell_centered_flag_is_read(case001_zip):
    smv = parse(case001_zip, 'case001')
    flags = [rec['CELL_CENTERED']
             for rec in smv['files']['SLICES'].values()]
    assert any(flags), "case001 writes cell-centered slices (SLCC records)"


def test_slice_records_carry_quantity_and_units(case001_zip):
    smv = parse(case001_zip, 'case001')
    for name, record in smv['files']['SLICES'].items():
        assert name.endswith('.sf')
        assert isinstance(record['QUANTITY'], str)
        assert len(record['QUANTITY']) > 0
        assert isinstance(record['UNITS'], str)


def test_blocked_cells_in_plane(case001_zip):
    blocks = fds.getBlockedCellsInPlane(case001_zip, 'case001', 3, 0.0)
    assert blocks.ndim == 2
    assert set(np.unique(blocks)).issubset({0.0, 1.0})
    assert blocks.sum() > 0, "the floor plane of case001 is obstructed"
