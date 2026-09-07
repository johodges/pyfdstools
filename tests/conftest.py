"""Shared fixtures for the pyfdstools test suite.

Every test runs against the FDS cases bundled in
``pyfdstools/examples/data``, so the suite needs no FDS installation and
no network access.
"""

import os
import shutil
import zipfile

import matplotlib
import pytest

# Select a non-interactive backend before pyplot is imported anywhere so
# that plotting tests do not require a display.
matplotlib.use('Agg')

DATA_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    'pyfdstools', 'examples', 'data')


def archive(chid):
    """Returns the path to the bundled archive for a case."""
    return os.path.join(DATA_DIR, '%s.zip' % (chid))


@pytest.fixture(scope='session')
def data_dir():
    """Directory holding the bundled example cases."""
    return DATA_DIR


@pytest.fixture(scope='session')
def case001_zip():
    """Archive of case001, a single mesh case with slice and boundary data."""
    return archive('case001')


@pytest.fixture(scope='session')
def case002_zip():
    """Archive of case002, a four mesh case with geometry output."""
    return archive('case002')


@pytest.fixture(scope='session')
def stretched_zip():
    """Archive of a case whose meshes use FDS mesh stretching."""
    return archive('stretched_mesh_example')


@pytest.fixture(scope='session')
def hfg_zip():
    """Archive of a four mesh case with heat flux gauge slices."""
    return archive('hfg_slice')


@pytest.fixture(scope='session')
def case001_dir(tmp_path_factory):
    """case001 extracted to a directory.

    Several readers take different code paths for a directory than for a
    zip archive (seek-based reads, glob-based file discovery), so both
    need coverage.
    """
    target = tmp_path_factory.mktemp('case001_dir')
    with zipfile.ZipFile(archive('case001'), 'r') as z:
        z.extractall(target)
    shutil.copy(os.path.join(DATA_DIR, 'case001.fds'),
                os.path.join(target, 'case001.fds'))
    return str(target)


@pytest.fixture()
def outdir(tmp_path):
    """Empty directory for files a test writes."""
    return str(tmp_path)
