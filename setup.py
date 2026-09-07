from setuptools import setup
import re

VERSIONFILE = "pyfdstools/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

with open("README.md", "rt", encoding="utf-8") as f:
    long_description = f.read()

# Optional third party libraries. These back the pyfdstools.paraview and
# pyfdstools.vtkhdf_plugin submodules, which are imported lazily so that
# the core package installs and imports without them.
paraview_requires = ['numpy-stl>=1.0', 'pyevtk>=1.0']
vtk_requires = ['vtk>=9.0']

setup(
    name='pyfdstools',
    version=verstr,
    description=(
        'This software is part of a python library to assist in developing '
        'and analyzing simulation results from Fire Dynamics Simulator '
        '(FDS). FDS is an open source computational fluid dynamics (CFD) '
        'software package developed by NIST. The sourcecode is available '
        'at: https://github.com/firemodels/fds'
    ),
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/johodges/pyfdstools',
    author='Jonathan Hodges',
    author_email='johodges@vt.edu',
    license='MIT',
    packages=['pyfdstools'],

    include_package_data=True,
    python_requires='>=3.8',
    install_requires=[
        'matplotlib>=3.0',
        'numpy>=1.17',
        'pandas>=0.25',
        'scipy>=1.3.1',
    ],
    extras_require={
        'paraview': paraview_requires,
        'vtk': vtk_requires,
        'all': paraview_requires + vtk_requires,
        'test': ['pytest>=7.0'],
    },

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Topic :: Scientific/Engineering',
    ],
)
