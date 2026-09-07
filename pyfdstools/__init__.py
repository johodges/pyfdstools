#-----------------------------------------------------------------------
# Copyright (C) 2020, All rights reserved
#
# Jonathan L. Hodges
#
#-----------------------------------------------------------------------
#=======================================================================
#
# DESCRIPTION:
# This software is part of a python library to assist in developing and
# analyzing simulation results from Fire Dynamics Simulator (FDS).
# FDS is an open source software package developed by NIST. The source
# code is available at: https://github.com/firemodels/fds
#
# EXAMPLES:
# See the examples subroutine for example operation.
#
# NOTE ON OPTIONAL DEPENDENCIES:
# The core of pyfdstools only requires numpy, scipy, matplotlib and
# pandas. Two submodules wrap third party visualization libraries which
# are large and are therefore optional:
#
#     pyfdstools.paraview       requires pyevtk and numpy-stl
#     pyfdstools.vtkhdf_plugin  requires vtk
#
# These submodules are imported lazily. Importing pyfdstools succeeds
# without them; the ImportError is only raised when one of the names
# they provide is actually requested. Install them with:
#
#     pip install pyfdstools[paraview]
#     pip install pyfdstools[vtk]
#     pip install pyfdstools[all]
#
#=======================================================================
# # IMPORTS
#=======================================================================

import importlib as _importlib

from .smokeviewParser import *
from .extractS3D import *
from .extractBoundaryData import *
from .extractGeomData import *
from .extractPlot3Ddata import *
from .extractParticleData import *
from .fdsErrorCalculation import *
from .fdsTypes import *
from .fdsFileOperations import *
from .extractCSVdata import *
from .inputFileVerification import *
from .utilities import *
from .colorSchemes import *
from .examples import *
from ._version import __version__

#=======================================================================
# # OPTIONAL SUBMODULES
#=======================================================================

# Maps each optional submodule to the pip extra which provides its
# third party requirements. Names defined by these submodules are
# resolved on first access by the module level __getattr__ below.
_OPTIONAL_MODULES = {
    'paraview': 'paraview',
    'vtkhdf_plugin': 'vtk',
}

# Cache of {public name: submodule name} populated the first time an
# optional submodule is successfully imported.
_optional_cache = {}

# Records the ImportError raised by each optional submodule which could
# not be loaded, so that a helpful message can be appended to the
# AttributeError raised for an unknown name.
_optional_errors = {}


def _load_optional_module(moduleName):
    """Imports an optional submodule and caches the names it defines.

    Parameters
    ----------
    moduleName : str
        Name of the submodule to import, relative to pyfdstools

    Returns
    -------
    module
        The imported submodule

    Raises
    ------
    ImportError
        If a third party requirement of the submodule is missing. The
        message names the pip extra which installs it.
    """

    try:
        module = _importlib.import_module('.%s' % (moduleName), __name__)
    except ImportError as err:
        _optional_errors[moduleName] = err
        extra = _OPTIONAL_MODULES[moduleName]
        raise ImportError(
            "pyfdstools.%s requires optional dependencies which are not "
            "installed (%s). Install them with "
            "'pip install pyfdstools[%s]'." % (moduleName, err, extra)
        ) from err
    for name in dir(module):
        if not name.startswith('_'):
            _optional_cache.setdefault(name, moduleName)
    globals()[moduleName] = module
    return module


def __getattr__(name):
    """Resolves names provided by the optional submodules on demand.

    This implements PEP 562 module level attribute access so that
    ``import pyfdstools`` never fails because vtk or pyevtk are absent,
    while ``pyfdstools.exportSl3dDataToVtk`` still works when they are
    installed.

    A submodule which cannot be imported is skipped rather than allowed
    to raise, so that an unknown attribute still produces AttributeError
    (which ``hasattr`` and the interactive completers rely on). The
    names of the skipped submodules are appended to the error message so
    that a missing optional dependency is still easy to diagnose.

    Parameters
    ----------
    name : str
        Attribute being looked up on the pyfdstools package

    Returns
    -------
    object
        The requested attribute

    Raises
    ------
    ImportError
        If the name is requested by explicit submodule name and that
        submodule's dependencies are missing
    AttributeError
        If no optional submodule provides the requested name
    """

    if name in _OPTIONAL_MODULES:
        return _load_optional_module(name)

    for moduleName in _OPTIONAL_MODULES:
        try:
            module = _load_optional_module(moduleName)
        except ImportError:
            continue
        if hasattr(module, name):
            value = getattr(module, name)
            globals()[name] = value
            return value

    message = "module %s has no attribute %s" % (__name__, name)
    if len(_optional_errors) > 0:
        skipped = ', '.join(sorted(_optional_errors.keys()))
        message = message + (
            ". Note the optional submodule(s) %s could not be imported; "
            "install the optional dependencies with "
            "'pip install pyfdstools[all]' if the name is provided by "
            "one of them." % (skipped))
    raise AttributeError(message)


def __dir__():
    """Lists the package attributes, including loadable optional names.

    Returns
    -------
    list
        Sorted list of attribute names available on pyfdstools
    """

    names = set(globals().keys())
    for moduleName in _OPTIONAL_MODULES:
        try:
            module = _load_optional_module(moduleName)
        except ImportError:
            continue
        names.update([n for n in dir(module) if not n.startswith('_')])
    return sorted(names)
