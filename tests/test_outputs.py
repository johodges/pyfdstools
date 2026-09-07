"""Tests for the remaining FDS output readers and the plotting helpers."""

import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pytest

import pyfdstools as fds


# ---------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------

def test_load_csv_devc(case001_zip):
    data = fds.load_csv(case001_zip, 'case001', '_devc')
    assert 'Time' in data.columns
    assert len(data) > 0
    assert np.all(np.diff(data['Time'].values) > 0)


def test_load_csv_hrr(case001_zip):
    data = fds.load_csv(case001_zip, 'case001', '_hrr')
    assert 'Time' in data.columns
    assert 'HRR' in data.columns
    assert data['HRR'].max() > 0


def test_load_csv_zip_and_dir_agree(case001_zip, case001_dir):
    fromZip = fds.load_csv(case001_zip, 'case001', '_devc')
    fromDir = fds.load_csv(case001_dir, 'case001', '_devc')
    assert list(fromZip.columns) == list(fromDir.columns)
    assert np.allclose(fromZip.values, fromDir.values)


def test_load_csv_missing_suffix_raises(case001_zip):
    with pytest.raises(FileNotFoundError, match='No csv file'):
        fds.load_csv(case001_zip, 'case001', '_no_such_output')


# ---------------------------------------------------------------------
# Plot3D output
# ---------------------------------------------------------------------

def test_readPlot3Ddata(case001_zip):
    grid, data = fds.readPlot3Ddata('case001', case001_zip, 30.0)
    assert grid.shape[3] == 3
    assert data.shape[3] == 5
    assert data.shape[:3] == grid.shape[:3]
    temperature = data[:, :, :, 0]
    finite = temperature[np.isfinite(temperature)]
    assert finite.size > 0
    assert finite.min() > -50 and finite.max() < 2000


def test_readPlot3Ddata_picks_the_nearest_time(case001_zip):
    """Different query times must select different plot3D files."""
    _, early = fds.readPlot3Ddata('case001', case001_zip, 10.0)
    _, late = fds.readPlot3Ddata('case001', case001_zip, 120.0)
    assert not np.allclose(early, late, equal_nan=True)


def test_readPlot3Ddata_missing_case_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        fds.readPlot3Ddata('nosuchcase', str(tmp_path), 0.0)


def test_findSliceLocation_from_plot3d(case001_zip):
    grid, data = fds.readPlot3Ddata('case001', case001_zip, 30.0)
    x, z, T, U, V, W, HRR = fds.findSliceLocation(
        grid, data, 2, 4.4, plot3d=True)
    assert x.shape == z.shape == T.shape == HRR.shape


# ---------------------------------------------------------------------
# Particle output
# ---------------------------------------------------------------------

def test_importParticle(case001_zip):
    prt5Files = fds.getFileList(case001_zip, 'case001', 'prt5')
    assert len(prt5Files) > 0
    result = fds.importParticle(prt5Files[0])
    assert result is not None


# ---------------------------------------------------------------------
# Smoke3D output
# ---------------------------------------------------------------------

def test_extractS3dValues(case001_zip):
    values, times = fds.extractS3dValues(case001_zip, 'case001')
    assert values is not None and len(values) > 0
    assert len(times) > 0
    for meshNum in values:
        for quantity, data in values[meshNum].items():
            assert isinstance(quantity, str)
            assert np.asarray(data).size > 0


def test_extractS3dValues_missing_files(tmp_path, capsys):
    values, times = fds.extractS3dValues(str(tmp_path), 'nosuchcase')
    assert values is None and times is None


# ---------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------

def test_plotSlice_returns_figure_and_axes(case001_zip):
    data, units = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    fig, ax = fds.plotSlice(
        data['x'], data['z'], data['datas'][:, :, -1], 1,
        clabel='TEMPERATURE (%s)' % (units), qnty_mn=0, qnty_mx=1000)
    assert ax.get_xlabel() and ax.get_ylabel()
    plt.close(fig)


@pytest.mark.parametrize('extend', ['both', 'below', 'above', 'neither',
                                    'min', 'max'])
def test_plotSlice_accepts_every_extend_spelling(case001_zip, extend):
    """'below' and 'above' are pyfdstools names for 'min' and 'max'.

    Passing them straight to matplotlib raised a ValueError.
    """

    data, _ = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    fig, ax = fds.plotSlice(
        data['x'], data['z'], data['datas'][:, :, -1], 1, extend=extend)
    plt.close(fig)


def test_plotSlice_rejects_unknown_extend(case001_zip):
    data, _ = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    with pytest.raises(ValueError, match='extend must be one of'):
        fds.plotSlice(data['x'], data['z'], data['datas'][:, :, -1], 1,
                      extend='sideways')


def test_plotSlice_does_not_modify_its_input(case001_zip):
    """The image path masks out-of-range values; it must copy first.

    Writing the mask into the caller's array left NaNs in data the
    caller still owned.
    """

    data, units = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    frame = data['datas'][:, :, -1]
    before = frame.copy()
    fig, ax = fds.plotSlice(
        data['x'], data['z'], frame, 1,
        contour=False, extend='below',
        qnty_mn=float(np.nanmin(frame)),
        qnty_mx=float(np.nanmedian(frame)))
    plt.close(fig)
    assert np.array_equal(before, frame, equal_nan=True)


def test_plotSlice_accepts_smv_colormap(case001_zip):
    data, _ = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    fig, ax = fds.plotSlice(
        data['x'], data['z'], data['datas'][:, :, -1], 1, cmap='SMV')
    plt.close(fig)


def test_renderSliceCsvs(case001_zip, outdir):
    data, units = fds.query2dAxisValue(
        case001_zip, 'case001', 'TEMPERATURE', 1, 2.55, time=30, dt=60)
    fds.renderSliceCsvs(data, 'case001', outdir)
    written = [f for f in os.listdir(outdir) if f.endswith('.csv')]
    assert len(written) == len(data['times'])
    import pandas as pd
    frame = pd.read_csv(os.path.join(outdir, written[0]), index_col=0)
    # Rows are the second in-plane coordinate, columns the first.
    assert frame.shape == data['datas'].shape[:2][::-1]
