"""Regression tests for bugs reported on the issue tracker.

Each test names the issue it covers and asserts on the behaviour that
was wrong, so that a reappearance is reported as that issue rather than
as an unexplained failure.
"""

import os
import shutil
import zipfile

import numpy as np
import pytest

import pyfdstools as fds
from pyfdstools.extractPlot3Ddata import _listSliceFilesForMesh


# =====================================================================
# Issue #8: readP3Dfile returns incorrect data
# https://github.com/johodges/pyfdstools/issues/8
# =====================================================================

def plot3dFile(workingDir, chid, timeStr):
    matches = [x for x in fds.getFileList(workingDir, chid, 'q')
               if timeStr in os.path.basename(x)]
    assert len(matches) == 1
    return matches[0]


def test_readP3Dfile_starts_at_the_data_record(case001_zip):
    """The values follow three header records, not the file start.

    readP3Dfile read from offset 0, so the first twelve values it
    returned were header bytes reinterpreted as floats and everything
    after them was shifted by twelve positions.
    """

    qfile = plot3dFile(case001_zip, 'case001', '120p00')

    f = fds.zopen(qfile)
    raw = f.read()
    f.close()

    nx, ny, nz = np.frombuffer(raw, dtype=np.int32, count=5)[1:4]
    count = int(nx)*int(ny)*int(nz)*5

    # The Fortran record marker before the values states their size, and
    # the file is exactly the header plus the values plus the closing
    # marker. Together these fix where the values begin.
    declared = int(np.frombuffer(raw, dtype=np.int32, count=1, offset=44)[0])
    assert declared == count*4
    assert len(raw) == 48 + count*4 + 4

    data, header = fds.readP3Dfile(qfile)
    assert list(header) == [nx, ny, nz]
    expected = np.frombuffer(raw, dtype=np.float32, count=count, offset=48)
    assert np.array_equal(data.flatten(order='F'), expected)


def test_readP3Dfile_returns_physical_values(case001_zip):
    """The decoded quantities must be physically possible.

    Reading from the wrong offset carried the tail of each quantity into
    the start of the next, which turned a 75 C temperature into a
    75 m/s velocity in this case.
    """

    data, header = fds.readP3Dfile(plot3dFile(case001_zip, 'case001', '120p00'))
    temperature, u, v, w, hrrpuv = [data[:, j] for j in range(0, 5)]

    # Ambient in this case is 20 C and nothing is cooled below it.
    assert temperature.min() >= 20.0 - 1e-3
    assert temperature.max() < 2000.0
    # A compartment fire of this size has velocities of a few m/s.
    for component in (u, v, w):
        assert np.abs(component).max() < 20.0
    # HRRPUV is a source term.
    assert hrrpuv.min() >= -1e-6


def test_readP3Dfile_rejects_a_malformed_file(case001_zip, outdir):
    f = fds.zopen(plot3dFile(case001_zip, 'case001', '120p00'))
    raw = f.read()
    f.close()

    # Truncated part way through the values.
    truncated = os.path.join(outdir, 'truncated.q')
    with open(truncated, 'wb') as f:
        f.write(raw[:len(raw)//2])
    with pytest.raises(ValueError, match='truncated'):
        fds.readP3Dfile(truncated)

    # Shorter than the header itself.
    stub = os.path.join(outdir, 'stub.q')
    with open(stub, 'wb') as f:
        f.write(raw[:20])
    with pytest.raises(ValueError, match='too short'):
        fds.readP3Dfile(stub)

    # A data record marker that disagrees with the grid header.
    corrupt = os.path.join(outdir, 'corrupt.q')
    with open(corrupt, 'wb') as f:
        f.write(raw[:44] + np.int32(12345).tobytes() + raw[48:])
    with pytest.raises(ValueError, match='does not look like a plot3D file'):
        fds.readP3Dfile(corrupt)


def test_writeP3Dfile_round_trips(outdir):
    rng = np.random.default_rng(0)
    reference = (rng.random((7, 5, 3, 5)) * 100).astype(np.float32)

    path = os.path.join(outdir, 'roundtrip.q')
    fds.writeP3Dfile(path, reference)

    data, header = fds.readP3Dfile(path)
    assert list(header) == [7, 5, 3]
    assert np.array_equal(
        data, np.reshape(reference, (7*5*3, 5), order='F'))


def test_writeP3Dfile_writes_valid_record_markers(outdir):
    """The file must carry the markers FDS and smokeview expect.

    Earlier releases wrote seven arbitrary floats where the second
    record's markers belong and omitted the data record's closing
    marker, leaving the file four bytes short.
    """

    data = np.zeros((4, 3, 2, 5), dtype=np.float32)
    path = os.path.join(outdir, 'markers.q')
    fds.writeP3Dfile(path, data)

    raw = open(path, 'rb').read()
    i32 = lambda o: int(np.frombuffer(raw, dtype=np.int32, count=1, offset=o)[0])

    assert i32(0) == 12 and i32(16) == 12          # grid record
    assert [i32(4), i32(8), i32(12)] == [4, 3, 2]
    assert i32(20) == 16 and i32(40) == 16         # four reals
    assert i32(44) == data.size*4                  # data record opens
    assert i32(48 + data.size*4) == data.size*4    # and closes
    assert len(raw) == 48 + data.size*4 + 4


def test_writeP3Dfile_rejects_wrong_shape(outdir):
    with pytest.raises(ValueError, match='must have shape'):
        fds.writeP3Dfile(os.path.join(outdir, 'bad.q'),
                         np.zeros((4, 3, 2, 4), dtype=np.float32))


# =====================================================================
# Issue #6: loading 2-D slices from a directory differs from a zip
# https://github.com/johodges/pyfdstools/issues/6
# =====================================================================

@pytest.fixture(scope='module')
def extracted(tmp_path_factory):
    """Every bundled case extracted next to its archive."""
    root = tmp_path_factory.mktemp('extracted')
    dataDir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        'pyfdstools', 'examples', 'data')
    out = {}
    for chid in ['case001', 'case002', 'hfg_slice', 'stretched_mesh_example']:
        target = os.path.join(str(root), chid)
        os.makedirs(target, exist_ok=True)
        with zipfile.ZipFile(os.path.join(dataDir, '%s.zip' % chid)) as z:
            z.extractall(target)
        loose = os.path.join(dataDir, '%s.fds' % chid)
        if os.path.exists(loose):
            shutil.copy(loose, target)
        out[chid] = target
    return out


def test_mesh_selection_does_not_match_other_meshes(
        stretched_zip, extracted):
    """Mesh 1 must not pick up meshes 10 and above.

    The archive branch tested whether the mesh prefix appeared anywhere
    in the path, so in this 14 mesh case the files of meshes 10 to 14
    also matched mesh 1.
    """

    chid = 'stretched_mesh_example'
    for meshNum in [1, 2, 10, 14]:
        prefix = '%s_%d_' % (chid, meshNum)
        fromZip = [os.path.basename(x)
                   for x in _listSliceFilesForMesh(stretched_zip, chid, prefix)]
        fromDir = [os.path.basename(x)
                   for x in _listSliceFilesForMesh(extracted[chid], chid, prefix)]

        assert fromZip == fromDir
        assert len(fromZip) == 1
        assert fromZip[0].startswith(prefix)


def test_readSLCFquantities_accepts_a_directory_without_a_separator(
        case001_zip, extracted):
    """A directory path need not end in a path separator.

    The glob interpolated the directory and the CHID with nothing
    between them, so a path not already ending in a separator matched
    no files at all and the case appeared to contain no slices.
    """

    fromZip = fds.readSLCFquantities('case001', case001_zip)
    plain = extracted['case001'].rstrip(os.sep)

    for workingDir in (plain, plain + os.sep):
        quantities, files, dims, meshes, centers, units = \
            fds.readSLCFquantities('case001', workingDir)
        assert len(files) == len(fromZip[1]) > 0
        assert sorted(quantities) == sorted(fromZip[0])


@pytest.mark.parametrize('chid', ['case001', 'case002'])
def test_query2dAxisValue_agrees_between_zip_and_directory(
        chid, data_dir, extracted):
    axis, value = (1, 2.55) if chid == 'case001' else (3, 7.2)

    fromZip, unitsZip = fds.query2dAxisValue(
        os.path.join(data_dir, '%s.zip' % chid), chid, 'TEMPERATURE',
        axis, value, atol=1e-2)
    fromDir, unitsDir = fds.query2dAxisValue(
        extracted[chid], chid, 'TEMPERATURE', axis, value, atol=1e-2)

    assert fromZip is not None and fromDir is not None
    assert unitsZip == unitsDir
    assert np.allclose(fromZip['x'], fromDir['x'], equal_nan=True)
    assert np.allclose(fromZip['datas'], fromDir['datas'], equal_nan=True)


@pytest.mark.parametrize('chid', ['case001', 'case002',
                                  'stretched_mesh_example'])
def test_readSLCF3Ddata_agrees_between_zip_and_directory(
        chid, data_dir, extracted):
    fromZip = fds.readSLCF3Ddata(
        chid, os.path.join(data_dir, '%s.zip' % chid), 'TEMPERATURE')
    fromDir = fds.readSLCF3Ddata(chid, extracted[chid], 'TEMPERATURE')

    assert fromZip[1] is not False
    assert np.allclose(fromZip[0], fromDir[0], equal_nan=True)
    assert np.allclose(fromZip[1], fromDir[1], equal_nan=True)
    assert fromZip[3] == fromDir[3]


def test_getFileList_is_sorted_and_order_matches_between_sources(
        case002_zip, extracted):
    """Mesh order decides the value on a shared face, so fix it.

    glob returns directory order and a zip returns member order; the two
    disagreed, which made case002 differ by 9.8 C at the corner where
    its four meshes meet depending on which the results were read from.
    """

    for ext in ('sf', 'xyz'):
        fromZip = [os.path.basename(x)
                   for x in fds.getFileList(case002_zip, 'case002', ext)]
        fromDir = [os.path.basename(x)
                   for x in fds.getFileList(extracted['case002'], 'case002', ext)]
        assert fromZip == sorted(fromZip)
        assert fromDir == sorted(fromDir)
        assert fromZip == fromDir
