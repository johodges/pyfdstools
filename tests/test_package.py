"""Tests for package level behaviour: imports, version, optional deps."""

import importlib
import sys

import pytest

import pyfdstools as fds


def test_version_is_a_string():
    assert isinstance(fds.__version__, str)
    assert len(fds.__version__.split('.')) >= 2


def test_core_import_needs_no_optional_dependencies(monkeypatch):
    """Importing pyfdstools must not require vtk, pyevtk or numpy-stl.

    These are large optional dependencies. Regressions here are easy to
    introduce by adding an unconditional import to a submodule, and they
    break every user who installed only the core requirements.
    """

    blocked = ('vtk', 'vtkmodules', 'evtk', 'stl', 'cv2')

    class Blocker:
        def find_module(self, name, path=None):
            if name.split('.')[0] in blocked:
                return self
            return None

        def load_module(self, name):
            raise ImportError('blocked for test: %s' % (name))

    for name in list(sys.modules):
        if name == 'pyfdstools' or name.startswith('pyfdstools.'):
            monkeypatch.delitem(sys.modules, name, raising=False)
        elif name.split('.')[0] in blocked:
            monkeypatch.delitem(sys.modules, name, raising=False)

    monkeypatch.setattr(sys, 'meta_path', [Blocker()] + sys.meta_path)
    module = importlib.import_module('pyfdstools')
    assert hasattr(module, 'query2dAxisValue')
    assert hasattr(module, 'parseSMVFile')


def test_unknown_attribute_raises_attribute_error():
    """A missing optional dependency must not turn into an ImportError.

    hasattr() and the interactive completers rely on AttributeError, so
    the lazy loader has to swallow the ImportError from an optional
    submodule it cannot load.
    """

    with pytest.raises(AttributeError):
        fds.thisNameDoesNotExistAnywhere
    assert not hasattr(fds, 'thisNameDoesNotExistAnywhere')


def test_public_names_are_exported():
    expected = [
        'query2dAxisValue', 'queryBndf', 'readSLCF3Ddata', 'readPlot3Ddata',
        'parseSMVFile', 'plotSlice', 'load_csv', 'fdsFileOperations',
        'getFileList', 'getEndianness', 'readSLCFquantities',
        'importBoundaryFile', 'buildSMVcolormap',
    ]
    missing = [name for name in expected if not hasattr(fds, name)]
    assert missing == []
