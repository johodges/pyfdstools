# Installation

pyfdstools needs Python 3.8 or newer.

## From PyPI

```bash
python -m pip install pyfdstools
```

## From source

```bash
git clone https://github.com/johodges/pyfdstools
python -m pip install ./pyfdstools
```

Use an editable install if you intend to modify the code or track the
repository with `git pull`:

```bash
git clone https://github.com/johodges/pyfdstools
python -m pip install -e ./pyfdstools
```

## In a virtual environment

Recommended, so that pyfdstools and its dependencies stay separate from
your system Python.

```bash
python -m venv /path/to/myenv
source /path/to/myenv/bin/activate        # Linux and macOS
/path/to/myenv/Scripts/activate           # Windows
python -m pip install pyfdstools
```

## Dependencies

The core package needs only four packages, all of which pip installs for
you:

| Package | Minimum | Used for |
| --- | --- | --- |
| numpy | 1.17 | every array operation |
| scipy | 1.3.1 | interpolation between meshes, convex hulls |
| matplotlib | 3.0 | plotting and colormaps |
| pandas | 0.25 | csv output as data frames |

Everything documented in this wiki works with just those, including
reading every FDS binary output format.

## Optional extras

Two submodules wrap larger third party libraries. They are loaded on
demand, so importing pyfdstools works whether or not they are installed;
the ImportError only appears when you reach for a name they provide, and
it names the extra to install.

| Extra | Installs | Enables |
| --- | --- | --- |
| `paraview` | `pyevtk`, `numpy-stl` | [Exporting to ParaView](Exporting-to-ParaView), `obstToStl` |
| `vtk` | `vtk` | reading FDS VTKHDF output (`pyfdstools.vtkhdf_plugin`) |
| `all` | both of the above | |
| `test` | `pytest` | running the test suite |

```bash
python -m pip install "pyfdstools[paraview]"
python -m pip install "pyfdstools[all]"
```

If you install pyfdstools and later find you want ParaView export,
installing the extra is enough; nothing in your existing scripts has to
change.

## Verifying the installation

```bash
python -c "import pyfdstools as fds; print(fds.__version__)"
```

Then run the bundled examples, which read the FDS cases shipped inside
the package and write their output to
`<install location>/pyfdstools/examples/generated`:

```bash
python -c "import pyfdstools as fds; fds.runExamples()"
```

`runExamples` returns a list of `(script, returncode, stdout, stderr)`
tuples, and takes `raiseOnError=True` if you want a failure to raise
rather than print.

## Blender

An older workflow drove Blender from pyfdstools; see
`pyfdstools/blender_examples.py` and the commented-out section of the
README. It requires Blender's bundled Python to be replaced with an
environment that has pyfdstools installed, and is not currently tested.
